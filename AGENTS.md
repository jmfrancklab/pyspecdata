# Repository Instructions

- Follow [STYLE_GUIDE.md](STYLE_GUIDE.md) for pyspecdata-native analysis,
  examples, and notebooks. Other projects and tools can reference that guide
  directly.
- Locally, before making changes or running tests, activate the shared base
  environment with `source /home/jmfranck/base/bin/activate`. Use only this
  environment; the checkout is already installed in editable mode. Do not
  reinstall unless explicitly asked.
- On cloud, before making changes or running tests, install from the repository
  root with `pip install -e . --no-build-isolation`. Install missing build,
  runtime, or testing dependencies as needed and retry if installation fails.
- Validate changes with the relevant `pytest` targets in the configured
  environment, and follow the linting and validation rules in `.github`.
