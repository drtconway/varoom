FROM ubuntu AS builder
ADD misc/docker/gather_shared.sh /gather_shared.sh
ADD include/ include/
ADD src/ src/
ADD tests/ tests/
ADD jamroot.jam jamroot.jam
RUN apt update && \
    apt upgrade -y && \
    apt install -y gcc g++ \
        libboost-dev libboost-iostreams-dev libboost-log-dev libboost-program-options-dev libboost-test-dev libboost-tools-dev \
        libhts-dev \
        nlohmann-json-dev && \
    b2 release && \
    cp ./bin/gcc-7/release/klbam /bin && \
    /gather_shared.sh /bin/klbam

FROM ubuntu
ADD misc/docker/unpack_shared.sh /unpack_shared.sh
COPY --from=builder /.shared_objects.tgz /
RUN /unpack_shared.sh && rm /unpack_shared.sh
