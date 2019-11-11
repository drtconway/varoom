FROM ubuntu
RUN apt update && \
    apt upgrade -y && \
    apt install -y gcc g++ git \
        libboost-dev libboost-iostreams-dev libboost-log-dev libboost-program-options-dev libboost-test-dev libboost-tools-dev \
        libhts-dev \
        nlohmann-json-dev && \
    git clone https://github.com/drtconway/varoom.git  && \
    cd varoom  && \
    b2 release  && \
    cp ./bin/gcc-7/release/klbam /bin  && \
    cd ..  && \
    rm -rf varoom && \
    apt remove -y gcc g++ git && \
    apt autoremove -y && \
    rm -rf /var/lib/apt/lists/*
