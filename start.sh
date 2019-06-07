#!/bin/sh

cp /app/data/config.yaml /app
snakemake -j 24
