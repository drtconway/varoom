#!/bin/bash

sofn=/.shared_object_list.txt
tmpfn=/.shared_object_list-tmp.txt
sotar=/.shared_objects.tgz

rm -f ${sofn}

for f in $@
do
    echo $f >> ${sofn}
    ldd ${f} | awk '{print $3}' >> ${sofn}
done

sort -u < ${sofn} | grep -v '^$' > ${tmpfn}
mv ${tmpfn} ${sofn}

tar -czvh -f ${sotar} -T ${sofn}

rm -f ${sofn}
