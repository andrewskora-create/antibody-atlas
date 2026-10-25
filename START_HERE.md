# Start Antibody Atlas

**Preview v0.1.0-rc4.** A local research workbench with guided setup, not a signed standalone application. You do not need a paid prediction account, Claude, ChatGPT or a GPU.

## Read the illustrated guide

**[Download the editable Word guide directly](https://github.com/andrewskora-create/antibody-atlas/releases/download/v0.1.0-rc4/Antibody_Atlas_User_Guide.docx)**. No installation is needed to read it.

After extracting the release ZIP, open **static/user-guide.html** in your browser. It works before the server is running. Once Atlas is running, choose **Quick start & user guide** at the top of the app. An editable Word copy is in **static/docs/Antibody_Atlas_User_Guide.docx**. The guide covers setup, the first example, substitutions, panels, exports, privacy and troubleshooting.

## Mac quick start

1. Download **antibody-atlas-0.1.0-rc4.zip** from [GitHub Releases](https://github.com/andrewskora-create/antibody-atlas/releases/tag/v0.1.0-rc4). Extract it and keep the whole folder together.
2. Install [Homebrew](https://brew.sh/) if needed, including its final shell setup instructions. In Terminal, run:

   ```sh
   brew install python@3.12 hmmer git
   ```

   If build tools are requested, run `xcode-select --install` and complete installation.
3. Double-click **Start Antibody Atlas.command** in the extracted folder. Use macOS's normal approval for this specific downloaded script if requested; do not disable Gatekeeper globally.
4. Leave Terminal open. First setup installs dependencies and models and needs internet. Allow several minutes or longer and several GB of free disk space.
5. When the browser opens, click **Load example**, then **Screen antibody**. Confirm the structure step finishes before interpreting structural results. Failed or unavailable evidence is not a pass.

Terminal alternative: type `cd `, drag the extracted folder into Terminal, press Return, then run `bash run_local.sh`. The default address is http://127.0.0.1:8000/.

## Linux and Windows

Use Ubuntu 24.04 with Python 3.11 or 3.12. Windows users need [WSL2](https://learn.microsoft.com/windows/wsl/install); run the commands inside Ubuntu, not PowerShell. Native Windows is not supported for the full structure workflow.

```sh
sudo apt update
sudo apt install python3 python3-venv python3-dev build-essential hmmer git
```

Enter the actual extracted Atlas folder, then run `bash run_local.sh`. Open the printed localhost address in Chrome. See [platform status](PLATFORMS.md) for validation limits. End-to-end Windows/WSL2 and Intel Mac operation have not been verified on those machines.

## Add clinical references after setup

In a second terminal in the extracted Atlas folder:

```sh
bash run_local.sh --setup-only --references
```

This explicitly builds the optional reference corpus from upstream data; see [data provenance](DATA_PROVENANCE.md). Save active results, stop and restart Atlas, then refresh. Until installed, clinical frequency is unavailable and substitution ranking uses residue similarity. PROPERMAB is a separate academic-only option, not part of default setup.

## Save and stop

Download results before stopping the server. Browser analyses are held in the session, not automatically archived. Use Ctrl+C in the server terminal to stop. Next time, use the same launcher. If port 8000 already has Atlas running, the launcher opens that session, which may be an older version. Save and stop it before switching versions, or run `PORT=8010 bash run_local.sh` for a separate local port.

## When something fails

Run `bash run_local.sh --doctor` from the extracted folder after setup. For an interrupted installation, check internet/disk space and rerun the launcher. For missing `hmmscan`, install HMMER. Read the detailed troubleshooting section in **static/user-guide.html** before deleting or replacing anything. Report issues with a public example and remove confidential data from logs and screenshots.

The core license permits academic and commercial use subject to its terms; optional tools and data retain their own terms. Atlas is provided as is, without warranty. Research predictions require experimental validation and do not establish manufacturing readiness.
