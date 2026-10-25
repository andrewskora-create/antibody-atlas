#!/usr/bin/env bash
# Build the bundled corresponding source; no network access. See START_HERE.md.
set -euo pipefail
cd "$(dirname "$0")/.."
root="$PWD"
if [[ ! -d native-source/cgal42 || ! -d native-source/nanoshaper ]]; then
  echo 'Native sources are in the shareable release ZIP. Extract that complete package first.' >&2; exit 2
fi
# Hints are portable; CMake finds distro Boost/GMP/MPFR or Homebrew's prefix.
extra=()
if command -v brew >/dev/null 2>&1; then extra+=("-DCMAKE_PREFIX_PATH=$(brew --prefix)"); fi
cmake -S native-source/cgal42 -B native-build/cgal -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DBoost_NO_BOOST_CMAKE=ON "${extra[@]}"
cmake --build native-build/cgal --parallel 2
cmake -S native-source/nanoshaper -B native-build/nanoshaper -DCGAL_DIR="$root/native-build/cgal" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DBoost_NO_BOOST_CMAKE=ON "${extra[@]}"
cmake --build native-build/nanoshaper --parallel 2
printf '\nBuild finished. Before starting the app, set:\nexport PROPERMAB_NANOSHAPER="%s/native-build/nanoshaper/NanoShaper"\n' "$root"
