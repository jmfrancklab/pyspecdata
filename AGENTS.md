# Repository Instructions

- When running locally: Before running tests or making code changes, always use the shared base environment at `/home/jmfranck/base` by running `source /home/jmfranck/base/bin/activate`. Do not reinstall the project locally unless explicitly asked; the local checkout is already installed in editable mode.
- When running on cloud: Before running tests or making code changes, always install the project in editable mode with Meson by running `pip install -e . --no-build-isolation` from the repository root.
- When running locally: Do not create or switch to any other virtual environment for this repository.
- Follow the linting and validation rules enforced by the workflows and scripts in `.github` when making changes.
- On cloud: Ensure all runtime and testing dependencies are installed so that `pytest` can execute without missing modules. If the editable install fails because build tools are absent (e.g., Meson, ninja, numpy), install the missing packages and rerun the command.
- On cloud: After installing dependencies, run the relevant `pytest` targets to validate changes.
- Locally: In the base environment, run the relevant `pytest` targets to validate changes.
