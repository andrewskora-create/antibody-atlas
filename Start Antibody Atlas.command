#!/usr/bin/env bash
# Double-click this file in Finder to start Antibody Atlas.
#
# macOS opens a .command file in Terminal and runs it. You do not have to type
# anything: this starts the local server and opens your browser. Leave the
# window open while you use the app, and close it (or press Ctrl-C) to stop.
#
# If Finder refuses to open it the first time, right-click the file and choose
# Open, then confirm. macOS blocks double-clicking a script it has marked as
# downloaded until you approve it once.

# Run from this file's own folder, whatever folder Finder was showing.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

printf '\033]0;Antibody Atlas\007'
echo "Antibody Atlas"
echo "Folder: $(pwd)"
echo

if [ ! -x ./run_local.sh ]; then
  if [ -f ./run_local.sh ]; then
    chmod +x ./run_local.sh
  else
    echo "ERROR: run_local.sh is not in this folder."
    echo "Keep this launcher inside the project folder."
    echo
    read -r -p "Press Return to close." _
    exit 1
  fi
fi

./run_local.sh
status=$?

echo
if [ "$status" -ne 0 ]; then
  echo "Antibody Atlas exited with an error (status $status)."
  echo "The messages above say what went wrong; scripts/doctor.py checks each engine."
else
  echo "Antibody Atlas has stopped."
fi
# Hold the window open so a failure is readable rather than vanishing.
read -r -p "Press Return to close this window." _
