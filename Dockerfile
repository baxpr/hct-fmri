# Start with FSL, ImageMagick, python3/pandas, xvfb base docker
FROM baxterprogers/fsl-base:v6.0.5.2

# Matlab reqs
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get -y install openjdk-8-jre && \
    apt-get clean
        
# Install the MCR
RUN wget -nv https://ssd.mathworks.com/supportfiles/downloads/R2023a/Release/6/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2023a_Update_6_glnxa64.zip \
    -O /opt/mcr_installer.zip && \
    unzip /opt/mcr_installer.zip -d /opt/mcr_installer && \
    /opt/mcr_installer/install -mode silent -agreeToLicense yes && \
    rm -r /opt/mcr_installer /opt/mcr_installer.zip

# Matlab env
ENV MATLAB_SHELL=/bin/bash
ENV AGREE_TO_MATLAB_RUNTIME_LICENSE=yes
ENV MATLAB_RUNTIME=/usr/local/MATLAB/MATLAB_Runtime/R2023a
ENV MCR_INHIBIT_CTF_LOCK=1
ENV MCR_CACHE_ROOT=/tmp

# Copy the pipeline code
COPY matlab /opt/hct-fmri/matlab
COPY src /opt/hct-fmri/src
COPY README.md /opt/hct-fmri/README.md

# Add pipeline to system path
ENV PATH /opt/hct-fmri/src:/opt/hct-fmri/matlab/bin:${PATH}

# Matlab executable must be run at build to extract the CTF archive
#RUN run_spm12.sh ${MATLAB_RUNTIME} function quit

# Entrypoint
ENTRYPOINT ["pipeline_entrypoint.sh"]
