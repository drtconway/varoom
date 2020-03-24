FROM drtomc/varoom-builder AS builder
ADD misc/docker/gather_shared.sh /gather_shared.sh
ADD include/ include/
ADD src/ src/
ADD tests/ tests/
ADD jamroot.jam jamroot.jam
RUN b2 release && \
    cp ./bin/gcc-7/release/klbam ./bin/gcc-7/release/varoom /bin && \
    /gather_shared.sh /bin/klbam /bin/varoom

FROM ubuntu
ADD misc/docker/unpack_shared.sh /unpack_shared.sh
COPY --from=builder /.shared_objects.tgz /
RUN /unpack_shared.sh && rm /unpack_shared.sh
