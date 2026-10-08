#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=10000M
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err
#SBATCH --job-name=sranges-validation
#SBATCH --array=0-199

SEED=$SLURM_ARRAY_TASK_ID 
cd "inf/${SEED}/"

JAR_PATH=$(realpath ../../ssr.jar)
VERSION_DIR=$(realpath ../../../Intellij_module_versions)

xvfb-run --auto-servernum  ~/skyline/zulu17.54.21-ca-fx-jdk17.0.13-linux_x64/bin/java -Djava.awt.headless=true -Djava.library.path=/lib/x86_64-linux-gnu/  \
    --add-modules javafx.fxml,javafx.base \
    -jar "$JAR_PATH" -threads -1 -seed 42 \
    -version_file "$VERSION_DIR/version-ssr.xml" \
    -version_file "$VERSION_DIR/version-sr.xml" \
    -version_file "$VERSION_DIR/version-feast.xml" \
    -version_file "$VERSION_DIR/version-beast.xml" \
    -version_file "$VERSION_DIR/version-sa.xml" \
    -version_file "$VERSION_DIR/version-beastLabs.xml" \
    -version_file "$VERSION_DIR/version-mm.xml" \
    -version_file "$VERSION_DIR/version-beastfx.xml" \
    -overwrite sRanges_inference.xml

cd ../../
