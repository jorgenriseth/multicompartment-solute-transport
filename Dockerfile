# Use Docker image with preinstalled conda
FROM continuumio/miniconda3
# Install ssh (missing dependency to run conda envs)
RUN apt-get update && \
    apt-get install -y ssh build-essential

# Upgrade conda
RUN conda upgrade -y conda

# Copy environment and requirements files into docker env
COPY environment.yml .
COPY requirements.txt .

# Update environment file with new environment name
RUN conda env update --file environment.yml --name dockerenv
SHELL ["conda", "run", "-n", "dockerenv", "/bin/bash", "-c"]

# Test dependencies
RUN python3 -c "import dolfin; print(dolfin.__version__); import h5py; print(h5py.__version__);import SVMTK; print(SVMTK)"

RUN echo "source activate dockerenv" > ~/.bashrc