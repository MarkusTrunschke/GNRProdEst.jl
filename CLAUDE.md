# Working in this repository

## Commits
- Do not add `Co-Authored-By` trailers to commit messages.

## Running things
- Bootstrap repetitions are threaded. Start Julia with `--threads=auto` when running
  anything that bootstraps, otherwise the repetitions run one after another.
- `boot_reps` defaults to 200 and every repetition re-estimates the whole model, so pass
  a small `boot_reps` when testing something out.

## Local files
- `run_example.jl` is a local scratch script and is gitignored. Do not commit it.
