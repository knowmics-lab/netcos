## Quick setup for working on this project with VS Code and GitHub

This repository contains the `drug_repurposing` codebase. Data files are large and are excluded from the git repo by default — use Git LFS or external storage for big datasets.

Follow the steps below (PowerShell snippets) to set up your local environment, Git, Git LFS, and push to GitHub.

1) Open the project in VS Code

- In VS Code: File -> Open Folder -> select the project folder (the folder that contains this README)
- Install recommended extensions if prompted (Python, Pylance, Jupyter, GitLens)

2) Create and activate the Conda environment (if you want to create fresh env)

```powershell
# create env from requirements (if you prefer conda env export, see notes below)
conda create -n drug_repurposing python=3.10 -y
conda activate drug_repurposing
pip install -r requirements.txt
```

If you prefer an `environment.yml` you can create one locally with:

```powershell
conda env export --name drug_repurposing --no-builds > environment.yml
```

3) Initialize git (if not already a git repo), commit and configure LFS

```powershell
cd "G:\My Drive\unict\unict 2024-25\drug_repurposing"
# if not already a repo:
git init
git add .gitattributes .gitignore
git commit -m "chore: add gitignore and gitattributes"

# Install Git LFS once on your machine (if not installed):
# https://git-lfs.github.com/ - or with scoop/choco on Windows
# Example (with git-lfs already in PATH):
git lfs install
# Track patterns (already present in .gitattributes, but ensure LFS is enabled):
git lfs track "data/**"
git lfs track "MITHrIL/input/**"
git add .gitattributes
git commit -m "chore: enable git lfs tracking for data"
```

4) Make a first commit of code (exclude big data as specified in .gitignore)

```powershell
git add src/ setup.py requirements.txt README-SETUP.md
git commit -m "feat: initial import of code and setup files"
```

5) Create a GitHub repo and push

Option A - Using HTTPS + PAT (recommended for new GitHub accounts):

```powershell
# create an empty repo on GitHub via website, then:
git remote add origin https://github.com/<your-username>/<repo-name>.git
# create a branch (e.g., main) and push
git branch -M main
git push -u origin main
```

Option B - Using SSH (if you already have an SSH key configured):

```powershell
git remote add origin git@github.com:<your-username>/<repo-name>.git
git push -u origin main
```

6) Using VS Code to work day-to-day

- Open the folder in VS Code.
- Use the Command Palette (Ctrl+Shift+P) -> Python: Select Interpreter -> choose the `drug_repurposing` conda env. If you activated the env before launching VS Code, `${env:CONDA_PREFIX}/python` setting will generally work.
- Use the Source Control view in VS Code to stage/commit changes and push/pull.

Notes and recommendations
- If datasets are too large for Git LFS or you prefer not to store them in the repo, host them elsewhere (cloud, institutional storage) and add download scripts to `data/` (e.g., `scripts/download_data.sh`).
- Consider adding `pre-commit` hooks for formatting and linting. A basic `pre-commit` config can run `black` and `ruff`.
- To reproduce environment reproducibly, create and commit an `environment.yml`:

```powershell
conda env export --name drug_repurposing --no-builds > environment.yml
git add environment.yml
git commit -m "chore: add conda environment file"
```

If you want, I can add a `pre-commit` template and a simple GitHub Actions workflow next.
