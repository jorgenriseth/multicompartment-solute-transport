# Use Docker image with preinstalled conda
FROM conda/miniconda3
# Install ssh (missing dependency to run conda envs)
RUN apt-get update && \
    apt-get install -y ssh
# Upgrade conda
RUN conda upgrade -y conda
# Copy environment file into docker env
COPY environment.yml .
# Update environment file with new environment name
RUN conda env update --file environment.yml --name dockerenv
SHELL ["conda", "run", "-n", "dockerenv", "/bin/bash", "-c"]
# Test dolfin
RUN python3 -c "import dolfin; print(dolfin.__version__); import h5py; print(h5py.__version);import SVMTK; print(SVMTK.__version__)"

