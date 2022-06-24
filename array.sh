#!/bin/bash

command=$(head -${SLURM_ARRAY_TASK_ID} CURRENT_DIR/commands.sh | tail -1)

srun $command
