#!/bin/bash
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

if [ ! -f CMakeLists.txt ]; then
  echo No CMakeLists.txt file in $PWD
  exit 1
fi

cmake . -DCMAKE_BUILD_TYPE=Release \
-DBUILD_SOLVER_AS_DLL=ON \
-DBUILD_CONTROL_AS_DLL=ON \
-DUSE_CONCURRENT_RECOVERY=ON \
-DUSE_SP_RECOVERY=ON
echo -e "\nGenerating source code documentation..."
make doc >& make_doc.log
if [ -s DoxyWarnings.log ]; then
  echo -e "\nHere are the warnings from doxygen:"
  cat DoxyWarnings.log
fi

# Check out previous version from the gh-pages branch and compare
# with the newly generated doc while ignoring generated time stamps.
# Since the index.html file also contains the build date, we need to
# exclude this file from the diff (assuming it is rarely changed).
# Therefore moving this file around.
git checkout gh-pages
mv docs/index.html index.old
mv doc/solver_html/index.html index.new
git rm -rf docs
mv doc/solver_html docs
mv index.old docs/index.html
git add docs/*.* docs/search/*.*
if git diff -I "^Generated on .* for FEDEM" HEAD --quiet; then
  echo No changes in source code documentation
  git reset --hard HEAD
  exit 1
else
  mv index.new docs/index.html
  git add docs/index.html
  exit 0
fi
