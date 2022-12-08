FROM quay.io/fenicsproject/stable:current
RUN apt-get update && apt-get upgrade -y && apt-get autoremove -y
COPY . /home/fenics/multirat
