# V1DD functional metrics

A Code Ocean pipeline that reads the V1DD NWB sessions, estimates per-ROI functional
metrics across stimulus families, and writes a versioned data asset.

Full documentation is in [`docs/`](docs/) — one page per analysis family, plus
[`docs/pipeline.md`](docs/pipeline.md) for how a run is structured and
[`docs/outputs.md`](docs/outputs.md) for the asset layout.

## Running it

In the capsule, "Reproducible Run" executes `code/run`. Before launching, record the
commit the asset is built from — a reproducible run copies `code/` without `.git`, so
the version cannot be derived there and the pipeline refuses to start without it:

```bash
git rev-parse HEAD > code/CODE_VERSION      # or set $V1DD_CODE_VERSION
```

To check the environment resolves without starting a multi-hour run:

```bash
code/run --check-env
```

## Configuration

| variable | default | meaning |
|---|---|---|
| `V1DD_CODE_VERSION` | — | commit stamped into the asset; overrides `code/CODE_VERSION` |
| `V1DD_INPUT_ASSET` | `/data/409828_V1DD_Filtered` | mounted NWB dataset |
| `V1DD_RESULTS_DIR` | `/results` | captured by Code Ocean |
| `V1DD_VALIDATION_DIR` | `/scratch/v1dd_metrics_validation` | discarded; checking the asset is not part of it |
| `V1DD_OUTPUT_TARGET` | `scratch` | set to `results` by the entry point, so interactive runs cannot fake an asset |

## Layout

```
code/run              capsule entry; puts code/src on PYTHONPATH
code/run_pipeline.py  version gate -> processing -> validation -> metadata
code/src/v1dd_metrics package: nwb, responses, families/, schema, provenance, metadata
code/validation       in-run integrity checks
code/tests            pytest
docs/                 pipeline and per-family documentation
```

## Local development

```bash
pip install -e code/src[dev]
python -m pytest
```

The NWB input is only mounted in Code Ocean, but a shipped asset is enough to exercise
most of the code offline — see [`docs/pipeline.md`](docs/pipeline.md).
