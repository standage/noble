FROM ubuntu:17.10
LABEL maintainer "Daniel Standage <daniel.standage@gmail.com>"

WORKDIR /work

# Install system packages
RUN apt-get update && apt-get install -y git build-essential libz-dev wget python3-dev python3-pip bwa samtools
RUN pip3 install --upgrade pip setuptools==32.0.0

# Install nuclmm
RUN pip3 install git+https://github.com/standage/nuclmm.git@bb80d1014cf3e2f69f27f8275f86e2cbd8e3fc12

# Install wgsim
RUN wget https://github.com/lh3/wgsim/archive/a12da3375ff3b51a5594d4b6fa35591173ecc229.tar.gz
RUN tar -xzf a12da3375ff3b51a5594d4b6fa35591173ecc229.tar.gz
RUN cd wgsim-a12da3375ff3b51a5594d4b6fa35591173ecc229/ && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm && cp wgsim /usr/local/bin/ && cd -
RUN rm -r a12da3375ff3b51a5594d4b6fa35591173ecc229.tar.gz wgsim-a12da3375ff3b51a5594d4b6fa35591173ecc229/

# Install lighter
RUN wget https://github.com/mourisl/Lighter/archive/df39031f8254f8351852f9f8b51b643475226ea0.tar.gz
RUN tar -xzf df39031f8254f8351852f9f8b51b643475226ea0.tar.gz
RUN cd Lighter-df39031f8254f8351852f9f8b51b643475226ea0/ && make && cp lighter /usr/local/bin/ && cd -
RUN rm -r df39031f8254f8351852f9f8b51b643475226ea0.tar.gz Lighter-df39031f8254f8351852f9f8b51b643475226ea0/

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/archive/v2.27.1.tar.gz
RUN tar -xzf v2.27.1.tar.gz
RUN cd bedtools2-2.27.1/ && make && cp bin/* /usr/local/bin/ && cd -
RUN rm -r v2.27.1.tar.gz bedtools2-2.27.1/


# Install kevlar
RUN pip3 install pysam==0.11.2.2 networkx==1.11 pandas==0.20.3 scipy==1.0.0
RUN pip3 install git+https://github.com/dib-lab/khmer.git@cd6dbdf21a6b1e4fc4806a205af2c5474d9996c2
RUN pip3 install git+https://github.com/dib-lab/kevlar.git@e6e9de07ca0acddfbba3ddc2e367b0496fbf6fd0

# Install SnakeMake for workflow execution
RUN pip3 install snakemake==4.4.0

CMD bash
