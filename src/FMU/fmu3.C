/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file fmu3.cpp
  \brief FMU implementation for FEDEM
  \details This file implements the subset of the functions declared in the file
  fmi3Functions.h from the FMI3 standard, necessary to run a FEDEM model as an
  FMU in a co-simulation context.
*/

#undef UNICODE

#include "fmi3Functions.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#if defined _MSC_VER
  #include <windows.h>
  typedef HINSTANCE LibHandle;
  std::string pathSeparator = "\\";
#elif defined __GNUC__
  #include <dlfcn.h>
  #include <cstring>
  typedef void* LibHandle;
  std::string pathSeparator = "/";
#endif


//Macro for printing debug info to terminal
#ifdef FMU_DEBUG
#define DEBUG_STDERR(x) { std::cerr << x << std::endl; }
#define DEBUG_STDOUT(x) { std::cout << x << std::endl; }
#else 
#define DEBUG_STDERR(x) {}
#define DEBUG_STDOUT(x) {}
#endif

typedef void (*DLPROC)();
typedef int (*DLPROC_INIT)(int,char**,const char*,
                           const double*,const int,
                           const double*,const int,
                           const double*,int*);
typedef int (*DLPROC_DONE)(bool);
typedef int (*DLPROC_GETSTATESIZE)();
typedef bool (*DLPROC_SETEXTFUNC)(int,double);
typedef bool (*DLPROC_SOLVENEXT)(int*);
typedef double (*DLPROC_EVALFUNC)(int,const char*,double,int*);
typedef bool (*DLPROC_SAVETRANSFORMATIONSTATE)(double*,const int);
typedef int (*DLPROC_RESTARTFROMSTATE)(const double*,const int,const int);


DLPROC getFuncAddress(const LibHandle& lib,
                      const std::string& fName)
{
#if defined _MSC_VER
  DLPROC address = (DLPROC)GetProcAddress(lib,fName.c_str());
#elif defined __GNUC__
  DLPROC address = (DLPROC)dlsym(lib,fName.c_str());
#else
#error "Platform not supported, neither _MSC_VER nor __GNUC__ defined"
#endif

  if (!address)
    DEBUG_STDERR(" *** Could not get function address for " + fName);
  
  return address;
}

DLPROC_INIT solverInit;
DLPROC_DONE solverDone;
DLPROC_GETSTATESIZE getStateSize;
DLPROC_GETSTATESIZE getTransformationStateSize;
DLPROC_SETEXTFUNC setExtFunc;
DLPROC_SOLVENEXT solveNext;
DLPROC_EVALFUNC evalFunc;
DLPROC_SAVETRANSFORMATIONSTATE saveTransformationState;
DLPROC_RESTARTFROMSTATE restartFromState;


std::string parentPath(std::string path)
{
  size_t pos = 0;
  while(path.find(pathSeparator,pos+1) != std::string::npos)
  {
    pos = path.find(pathSeparator,pos+1);
  }
  return path.substr(0,pos);
}

#ifdef __cplusplus
extern "C" {
#endif
  
  enum fmuStateCode
  {
    FMUSTART= 1u << 0,
    FMUEND = 1u << 1,
    FMUINSTANTIATED = 1u << 2,
    FMUINITIALIZATION = 1u << 3,
    FMUSTEPCOMPLETE= 1u << 4,
    FMUSTEPINPROGRESS = 1u << 5,
    FMUSTEPFAILED = 1u << 6,
    FMUSTEPCANCELED = 1u << 7,
    FMUTERMINATED = 1u << 8,
    FMUERROR = 1u << 9,
    FMUFATAL = 1u << 10
  };
  
  //State of internal fmu parameters and variables, as well as solver state for restart.
  struct modelState
  {
    std::string modelIdentifier;
    std::string modelGuid;
    int numReals;
    int numInputs;
    int numOutputs;
    int numParams;
    int numTransforms;
    
    int* fedemInputIndices;
    int* fedemOutputIndices;
    int* fedemTransformIndices;
    
    int solverStateSize;
    double* solverState;
    int transformationStateSize;
    double* transformationState;
    double t;
    double* reals;
  };
  
  struct componentInstance
  {
    fmuStateCode stateCode;
    modelState state;
    modelState initialState;
    fmi3String instanceName;
    fmi3Boolean logging = false;
    const fmi3CallbackFunctions *functions;
  };
  
  
  void readConfig(componentInstance* comp, std::string path)
  {
    DEBUG_STDOUT("Configpath: " + path);
    // TODO(RunarHR): readConfig: Make this more robust
    std::string line;
    std::ifstream confFile(path);
    std::getline(confFile,line);
    comp->initialState.modelIdentifier = line;
    std::getline(confFile,line);
    comp->initialState.modelGuid = line;
    std::getline(confFile,line);
    comp->initialState.numReals = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numInputs = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numOutputs = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numParams = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numTransforms = atoi(line.c_str());
    comp->initialState.fedemInputIndices = (int*)comp->functions->allocateMemory( comp->initialState.numInputs, sizeof(fmi3Integer) );
    comp->initialState.fedemOutputIndices = (int*)comp->functions->allocateMemory( comp->initialState.numOutputs, sizeof(fmi3Integer) );
    comp->initialState.fedemTransformIndices = (int*)comp->functions->allocateMemory( comp->initialState.numTransforms, sizeof(fmi3Integer) );
    comp->state.fedemInputIndices = (int*)comp->functions->allocateMemory( comp->initialState.numInputs, sizeof(fmi3Integer) );
    comp->state.fedemOutputIndices = (int*)comp->functions->allocateMemory( comp->initialState.numOutputs, sizeof(fmi3Integer) );
    comp->state.fedemTransformIndices = (int*)comp->functions->allocateMemory( comp->initialState.numTransforms, sizeof(fmi3Integer) );
    
    for(int i=0; i< comp->initialState.numInputs; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemInputIndices[i] = atoi(line.c_str());
    }
    for(int i=0; i< comp->initialState.numOutputs; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemOutputIndices[i] = atoi(line.c_str());
    }
    for(int i=0; i< comp->initialState.numTransforms; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemTransformIndices[i] = atoi(line.c_str());
    }
    
    confFile.close();
  }
  
  void copyModelState(modelState* destination, modelState* source)
  {
    //Copy arrays
    memcpy( destination->reals, source->reals, sizeof(double)*source->numReals );
    memcpy( destination->solverState, source->solverState, sizeof(double)*(source->solverStateSize ));
    memcpy( destination->transformationState, source->transformationState, sizeof(double)*(source->transformationStateSize ));
    
    //Copy variables
    destination->modelIdentifier = source->modelIdentifier;
    destination->modelGuid = source->modelGuid;
    destination->numReals = source->numReals;
    destination->numInputs = source->numInputs;
    destination->numOutputs = source->numOutputs;
    destination->numParams = source->numParams;
    destination->numTransforms = source->numTransforms;
    
    destination->t = source->t;
    destination->solverStateSize = source->solverStateSize;
    destination->transformationStateSize = source->transformationStateSize;
    
    memcpy( destination->fedemInputIndices, source->fedemInputIndices, sizeof(int)*source->numInputs );
    memcpy( destination->fedemOutputIndices, source->fedemOutputIndices, sizeof(int)*source->numOutputs );
    memcpy( destination->fedemTransformIndices, source->fedemTransformIndices, sizeof(int)*source->numTransforms );
  }
  
  fmi3Component fmi3Instantiate(fmi3String  instanceName, fmi3Type fmuType, fmi3String  GUID, fmi3String  fmuResourceLocation, const fmi3CallbackFunctions* functions, fmi3Boolean visible, fmi3Boolean loggingOn)
  {
    componentInstance* comp = 0;
    
    comp = (componentInstance*)functions->allocateMemory(1,sizeof(componentInstance));
    comp->logging = loggingOn;
    
    comp->functions = functions;
    comp->instanceName = instanceName;
    
    
    //Get resource folder path
#if defined _MSC_VER
    char path[1024];
    HMODULE hm = NULL;
    
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | 
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR) &fmi3GetTypesPlatform, 
                            &hm))
    {
      int ret = GetLastError();
      if(comp->logging)comp->functions->logger(functions->componentEnvironment, instanceName, fmi3Error, "error",
                                               "fmi3Instantiate: Could not retrieve path of running module from GetModuleHandleEx.");
    }
    GetModuleFileNameA(hm, path, sizeof(path));
#elif defined __GNUC__
    Dl_info info;
    if (dladdr((void*)fmi3Instantiate, &info))
    {
      DEBUG_STDOUT("FMU-library location: " + std::string(info.dli_fname));
    }
    else
    {
      DEBUG_STDERR(" *** Could not find location of FMU-library");
      return 0;
    }
    const char* path = info.dli_fname;
#endif
    
    std::string platformPath(parentPath(path));
    DEBUG_STDOUT("platformPath: " + platformPath);
    std::string binariesPath(parentPath(platformPath.c_str()));
    DEBUG_STDOUT("binariesPath: " + binariesPath);
    std::string rootPath(parentPath(binariesPath.c_str()));
    std::string resourcePath = rootPath + pathSeparator + "resources";
    DEBUG_STDOUT("resourcePath: " + resourcePath);

    readConfig(comp, resourcePath + pathSeparator + "config.txt");

    if(strcmp(comp->initialState.modelGuid.c_str(), GUID ) != 0)
    {
      DEBUG_STDOUT("fmi3Instantiate: GUIDs does not match");
      if(comp->logging)functions->logger(functions->componentEnvironment, instanceName, fmi3Error, "error",
                                         "fmi3Instantiate: GUIDs does not match.");
      DEBUG_STDOUT(comp->initialState.modelGuid.c_str());
      functions->freeMemory(comp);
      return 0;
    }

    std::string solverPath(getenv("FEDEM_SOLVER"));
    DEBUG_STDOUT("solverPath: " + solverPath);

    // NOTE(RunarHR): Initialization of solver should be done in EnterInitializationMode, but must be done here to get solverState size. 
#if defined _MSC_VER
    std::string solverBinariesPath(parentPath(solverPath.c_str()));
    SetDllDirectory(solverBinariesPath.c_str());
    LibHandle h_solver = LoadLibrary(solverPath.c_str());
    SetDllDirectory(NULL);
#elif defined __GNUC__
    LibHandle h_solver = dlopen(solverPath.c_str(), RTLD_LAZY);
#endif
    if (!h_solver) {
      DEBUG_STDERR(" *** Could not load solver library");
      return 0;
    }

    solverInit = (DLPROC_INIT)getFuncAddress(h_solver,"solverInit");
    solverDone = (DLPROC_DONE)getFuncAddress(h_solver,"solverDone");
    getStateSize = (DLPROC_GETSTATESIZE)getFuncAddress(h_solver,"getStateSize");
    getTransformationStateSize = (DLPROC_GETSTATESIZE)getFuncAddress(h_solver,"getTransformationStateSize");
    setExtFunc = (DLPROC_SETEXTFUNC)getFuncAddress(h_solver,"setExtFunc");
    solveNext = (DLPROC_SOLVENEXT)getFuncAddress(h_solver,"solveNext");
    evalFunc = (DLPROC_EVALFUNC)getFuncAddress(h_solver,"evalFunc");
    saveTransformationState = (DLPROC_SAVETRANSFORMATIONSTATE)getFuncAddress(h_solver,"saveTransformationState");
    restartFromState = (DLPROC_RESTARTFROMSTATE)getFuncAddress(h_solver,"restartFromState");

    std::string workingDir = resourcePath + pathSeparator + "model";

    // Lambda function adding option files to argvStart based on existance
    std::vector<char*> argvStart;
    auto&& addOptionFile=[workingDir,&argvStart](const char* opt,const char* fileName)
    {
      std::string filePath = workingDir + pathSeparator + fileName;
      std::ifstream fs(filePath.c_str());
      if (fs.good())
      {
        argvStart.push_back(const_cast<char*>(opt));
        argvStart.push_back(const_cast<char*>(fileName));
      }
    };

    argvStart.reserve(9);
    argvStart.push_back(const_cast<char*>("fedem_solver"));
    argvStart.push_back(const_cast<char*>("-cwd"));
    argvStart.push_back(const_cast<char*>(workingDir.c_str()));
    addOptionFile("-fco","fedem_solver.fco");
    addOptionFile("-fop","fedem_solver.fop");
    addOptionFile("-fao","fedem_solver.fao");

    DEBUG_STDOUT("Initializing solver");
    int status = solverInit(argvStart.size(),argvStart.data(),NULL,NULL,0,NULL,0,NULL,NULL);
    if (status < 0)
    {
      if(comp->logging) functions->logger(functions->componentEnvironment, instanceName, fmi3Error, "error",
                                          "Could not initialize solver.");
      DEBUG_STDERR(" *** Solver failed to initialize" + std::to_string(status));
      return 0;
    }
    
    comp->initialState.solverStateSize = getStateSize();
    comp->initialState.transformationStateSize = getTransformationStateSize();
    
    //Initial State
    comp->initialState.reals = (double*)functions->allocateMemory( comp->initialState.numReals, sizeof(fmi3Real) );
    comp->initialState.solverState = (double*)comp->functions->allocateMemory( comp->initialState.solverStateSize,sizeof(fmi3Real) );
    comp->initialState.transformationState = (double*)comp->functions->allocateMemory( comp->initialState.transformationStateSize,sizeof(fmi3Real) );
    comp->initialState.t = 0;
    
    //Working State
    comp->state.reals = (double*)functions->allocateMemory( comp->initialState.numReals, sizeof(fmi3Real) );
    comp->state.solverState = (double*)comp->functions->allocateMemory( comp->initialState.solverStateSize,sizeof(fmi3Real) );
    comp->state.transformationState = (double*)comp->functions->allocateMemory( comp->initialState.transformationStateSize,sizeof(fmi3Real) );
    comp->state.t = 0;
    
    copyModelState( &(comp->state), &(comp->initialState) );
    comp->stateCode = fmuStateCode::FMUINSTANTIATED;
    
    return comp;
  }
  
  //Stop simulation
  fmi3Status fmi3Terminate(fmi3Component c)
  {
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED))
    {
      //Close solver
      //solverDone(true);
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3DoStep(fmi3Component c, double currentCommunicationPoint, fmi3Real communicationStepSize, fmi3Boolean newStep)
  {
    DEBUG_STDOUT("fmi3DoStep");

    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & fmuStateCode::FMUSTEPCOMPLETE)
    {
      comp->state.t = currentCommunicationPoint;

      //SET INPUT FUNCTIONS IN FEDEM.
      int err = 0;
      for (int i = 0; i < comp->state.numInputs && err == 0; i++)
        if (!setExtFunc(comp->state.fedemInputIndices[i], comp->state.reals[i]))
          --err;

      if (err == 0) solveNext(&err);

      if (err != 0)
      {
        if(comp->logging)comp->functions->logger(comp->functions->componentEnvironment, comp->instanceName, fmi3Error, "error",
                                                 "fmi3DoStep: Solver step failed.");
        comp->stateCode = fmuStateCode::FMUSTEPFAILED;
        return fmi3Error;
      }

      //GET OUTPUT FROM FEDEM
      for (int i = 0; i < comp->state.numOutputs && err == 0; i++)
        comp->state.reals[comp->state.numInputs+i] = evalFunc(comp->state.fedemOutputIndices[i], NULL, -1.0, &err);

      if (err != 0)
      {
        comp->stateCode = fmuStateCode::FMUSTEPFAILED;
        return fmi3Error;
      }     

      //GET TRANSFORMATION OUTPUT FROM FEDEM.
      if(comp->state.numTransforms > 0)
      {
        saveTransformationState(comp->state.transformationState, comp->state.transformationStateSize);

        size_t idat = 3; //advance past timestep data
        //Loop over all objects in transformationState. Size per object is 14.
        double* mem = comp->state.transformationState;
        for (int i = 0; i < (comp->state.transformationStateSize-3)/14; i++)
        {
          idat++; //Skip typeId
          int baseID = (fmi3Integer)mem[idat++];

          //Check if baseID of transformation object match any of the
          //ones specified for the FMU.
          for (int j = 0; j < comp->state.numTransforms; j++)
            if (comp->state.fedemTransformIndices[j] == baseID)
            {
              int offset = comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + j*12;
              memcpy(comp->state.reals+offset, mem+idat, sizeof(double)*12);
              break;
            }

          idat += 12;
        }
      }

      comp->stateCode = fmuStateCode::FMUSTEPCOMPLETE;
      return fmi3OK;
    }
    
    // TODO(RunarHR): fmi3DoStep: If communicationStepSize is not equal to predefined step-size: return fmi3Error and write log-message
    
    comp->stateCode = fmuStateCode::FMUSTEPFAILED;
    return fmi3Error;
  }
  
  
  void fmi3FreeInstance(fmi3Component c)
  {
    //Clean up componentInstance, free all memory and resources
    
    if(c == 0)
      return;
    
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      //Free state
      comp->functions->freeMemory( comp->state.reals );
      comp->functions->freeMemory( comp->state.solverState );
      comp->functions->freeMemory( comp->state.transformationState );
      comp->functions->freeMemory( comp->state.fedemInputIndices );
      comp->functions->freeMemory( comp->state.fedemOutputIndices );
      
      //Free initial state
      comp->functions->freeMemory( comp->initialState.reals );
      comp->functions->freeMemory( comp->initialState.solverState );
      comp->functions->freeMemory( comp->initialState.transformationState );
      comp->functions->freeMemory( comp->initialState.fedemInputIndices );
      comp->functions->freeMemory( comp->initialState.fedemOutputIndices );
      
      comp->functions->freeMemory( comp );
      
      //Close solver
      solverDone(true);
    }
    
    return;
  }
  
  fmi3Status fmi3SetupExperiment(fmi3Component c, fmi3Boolean toleranceDefined, double tolerance, fmi3Real startTime, fmi3Boolean stopTimeDefined, fmi3Real stopTime)
  {
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode == fmuStateCode::FMUINSTANTIATED)
    {
      // TODO(RunarHR): fmi3SetupExperiment: setup start time etc. Called before initialize
      //Setup solver parameters
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3EnterInitializationMode(fmi3Component c) {
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode == fmuStateCode::FMUINSTANTIATED)
    {
      comp->stateCode = fmuStateCode::FMUINITIALIZATION;
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3ExitInitializationMode(fmi3Component c) {
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode == fmuStateCode::FMUINITIALIZATION)
    {
      // TODO(RunarHR): fmi3ExitInitializationMode: If values corresponding to input channels/external functions have been set, update solver. 
      
      comp->stateCode = fmuStateCode::FMUSTEPCOMPLETE;
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  
  //WARNING: This fails with Ansys TwinBuilder. The fmi3FMUstate pointer is not NULL on first call.
  fmi3Status fmi3GetFMUstate (fmi3Component c, fmi3FMUstate* FMUstate) {
    
    /*From Doc: fmi3GetFMUstate makes a copy of the internal FMU state and returns a pointer to this copy
    (FMUstate). If on entry *FMUstate == NULL, a new allocation is required. If *FMUstate !=
    NULL, then *FMUstate points to a previously returned FMUstate that has not been modified
    since. In particular, fmi3FreeFMUstate had not been called with this FMUstate as an argument.
    [Function fmi3GetFMUstate typically reuses the memory of this FMUstate in this case and
    returns the same pointer to it, but with the actual FMUstate.]*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state == 0)
      {
        state = (modelState*)comp->functions->allocateMemory(1,sizeof(modelState));
        state->reals = (double*)comp->functions->allocateMemory(comp->state.numReals,sizeof(fmi3Real));
        state->solverState = (double*)comp->functions->allocateMemory(comp->state.solverStateSize,sizeof(fmi3Real));
        state->transformationState = (double*)comp->functions->allocateMemory(comp->state.transformationStateSize,sizeof(fmi3Real));
      }
      
      copyModelState(state, &(comp->state));
      
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
    
  }
  fmi3Status fmi3SetFMUstate (fmi3Component c, fmi3FMUstate FMUstate) {
    /*From Doc: fmi3SetFMUstate copies the content of the previously copied FMUstate back and uses it as
    actual new FMU state. The FMUstate copy does still exist.*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state == 0)
      {
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi3Error;
      }
      
      //Copy input state to stored state
      copyModelState(&(comp->state), state);
      
      //Reset solver.
      restartFromState(comp->state.solverState, comp->state.solverStateSize,0); // NOTE(RunarHR): writeToRDB is 0. No results saved.
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  fmi3Status fmi3FreeFMUstate(fmi3Component c, fmi3FMUstate* FMUstate) {
    /* From Doc: fmi3FreeFMUstate frees all memory and other resources allocated with the fmi3GetFMUstate
    call for this FMUstate. The input argument to this function is the FMUstate to be freed. If a null
    pointer is provided, the call is ignored. The function returns a null pointer in argument FMUstate.
    These functions are only supported by the FMU, if the optional capability flag
    <fmiModelDescription> <ModelExchange / CoSimulation canGetAndSetFMUstate in =
    "true"> in the XML file is explicitly set to true (see sections 3.3.1 and 4.3.1).*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state != 0)
      {
        comp->functions->freeMemory(state->solverState);
        comp->functions->freeMemory(state->transformationState);
        comp->functions->freeMemory(state->reals);
        comp->functions->freeMemory(state);
      }
      
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3SetReal(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, const double value[])
  {
    DEBUG_STDOUT("fmi3SetReal");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE))
    {
      if(comp->state.numReals == 0)
      {
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi3Error;
      }

      for(size_t i = 0; i < nvr; i++)
      {
        if((int)vr[i] >= ( comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + comp->state.numTransforms))
        {
          comp->stateCode = fmuStateCode::FMUERROR;
          return fmi3Error;
        }
        comp->state.reals[vr[i]] = value[i];
      }
      
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetReal(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, double value[])
  {
    DEBUG_STDOUT("fmi3GetReal");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      DEBUG_STDOUT("stateCode-correct");
      if(comp->state.numReals == 0)
      {
        DEBUG_STDOUT("num reals == 0");
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi3Error;
      }

      for(size_t i = 0; i < nvr; i++)
      {
        if((int)vr[i] >= ( comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + comp->state.numTransforms*12 ))
        {
          DEBUG_STDOUT("Vr > numInputs + numOutputs + numParams + numTransforms*12");
          comp->stateCode = fmuStateCode::FMUERROR;
          return fmi3Error;
        }
        value[i] = comp->state.reals[vr[i]];
      }
      
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  const char* fmi3GetTypesPlatform() {
    DEBUG_STDOUT("fmi3GetTypesPlatform");
    return "default";
  }
  
  fmi3Status fmi3Reset(fmi3Component c)
  {
    DEBUG_STDOUT("fmi3Reset");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      //Copy initialState to state
      memcpy( comp->state.reals, comp->initialState.reals, sizeof(double)*(comp->initialState.numReals));
      memcpy( comp->state.solverState, comp->initialState.solverState, sizeof(double)*(comp->initialState.solverStateSize));
      memcpy( comp->state.transformationState, comp->initialState.transformationState, sizeof(double)*(comp->initialState.transformationStateSize));
      
      //Reset variables
      comp->state.t = comp->initialState.t;
      
      //Reset solver.
      //restartFromState(comp->state.solverState, comp->state.solverStateSize,0); // NOTE(RunarHR): writeToRDB is 0. No results saved.
      
      comp->stateCode = fmuStateCode::FMUINSTANTIATED;
      
      return fmi3OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  const char* fmi3GetVersion() {
    return fmi3Version;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //////////////////////////////////////////////////////
  //        NOT IMPLEMENTED                           //
  //////////////////////////////////////////////////////
  
  fmi3Status fmi3GetInteger(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, int value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetInteger
    DEBUG_STDOUT("fmiGetInteger");
    return fmi3OK;
  }
  
  
  fmi3Status fmi3GetBoolean(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, fmi3Boolean value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetBoolean
    DEBUG_STDOUT("fmiGetBoolean");
    return fmi3OK;
  }
  
  
  fmi3Status fmi3GetString(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, fmi3String value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetString
    return fmi3OK;
  }
  
  // Derivatives is not supported. 
  //canInterpolateInputs and MaxOutputDerivativeOrder must be set to false and 0 in element "CoSimulation" in the xml-file
  fmi3Status fmi3SetRealInputDerivatives(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, const int order[], const double value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3SetRealInputDerivatives
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  fmi3Status fmi3GetRealOutputDerivatives(fmi3Component c, const fmi3ValueReference vr[], size_t  nvr, const int order[], double value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetRealOutputDerivatives
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3CancelStep(fmi3Component c)
  {
    // TODO(RunarHR): [Optional] Implement fmi3CancelStep
    // Only relevant if doStep runs asynchronously.
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetStatus(fmi3Component c, const fmi3StatusKind s, fmi3Status* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetRealStatus(fmi3Component c, const fmi3StatusKind s, double* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetRealStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetIntegerStatus(fmi3Component c, const fmi3StatusKind s, int* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetIntegerStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetBooleanStatus(fmi3Component c, const fmi3StatusKind s, fmi3Boolean* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetBooleanStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3GetStringStatus(fmi3Component c, const fmi3StatusKind s, fmi3String* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi3GetStringStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3SetInteger(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, const int value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3SetInteger
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3SetBoolean(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, const fmi3Boolean value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3SetBoolean
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3SetString(fmi3Component c, const fmi3ValueReference vr[], size_t nvr, const fmi3String value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi3SetString
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
  fmi3Status fmi3SetDebugLogging(fmi3Component c, fmi3Boolean loggingOn, size_t nCategories, const fmi3String categories[]) 
  {
    // TODO(RunarHR): [Optional] Implement fmi3SetDebugLogging
    /*componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;*/
    return fmi3OK;
  }
  
  fmi3Status fmi3SerializedFMUstateSize(fmi3Component c, fmi3FMUstate FMUstate, size_t *size) {
    // TODO(RunarHR): [Optional] Implement fmi3SerializedFMUstateSize
    
    /*From Doc: fmi3SerializedFMUstateSize returns the size of the byte vector, in order that FMUstate can
    be stored in it. With this information, the environment has to allocate an fmi3Byte vector of the
    required length size.*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  fmi3Status fmi3SerializeFMUstate (fmi3Component c, fmi3FMUstate FMUstate, fmi3Byte serializedState[], size_t size) {
    // TODO(RunarHR): [Optional] Implement fmi3SerializeFMUstate
    
    /*From Doc: fmi3SerializeFMUstate serializes the data which is referenced by pointer FMUstate and
    copies this data in to the byte vector serializedState of length size, that must be provided by
    the environment.*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  fmi3Status fmi3DeSerializeFMUstate (fmi3Component c, const fmi3Byte serializedState[], size_t size,
                                      fmi3FMUstate* FMUstate) {
    // TODO(RunarHR): [Optional] Implement fmi3DeSerializeFMUstate
    
    /*From Doc: fmi3DeSerializeFMUstate deserializes the byte vector serializedState of length size,
    constructs a copy of the FMU state and returns FMUstate, the pointer to this copy. [The
    simulation is restarted at this state, when calling fmi3SetFMUState with FMUstate.]
    These functions are only supported by the FMU, if the optional capability flags
    canGetAndSetFMUstate and canSerializeFMUstate in
    <fmiModelDescription><ModelExchange / CoSimulation> in the XML file are explicitly set
    to true (see sections 3.3.1 and 4.3.1).*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi3Error;
  }
  
#ifdef __cplusplus
} // closing brace for extern "C"
#endif
