#include "varoom/sam.hpp"

#include <iostream>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sam_pileup tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    const char* sam_text =
		"NS500817:604:HLY3VBGXC:4:12403:17489:16364	99	18	29769769	60	75M	=	29769975	281	TAGAGAATGTATTCATTAAGAAGTCTCAGCCAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAGCCAACA	=:=?A>@>D>>?@DA?@=AEAAE?EAEAEE.BEDE6DAA?B@AE?AEADBDBEDABBBBBBBBBBB-AED.;<1=	MC:Z:75M	BD:Z:LLNLOMOLNJMLMLNLLKMLKMLONKLMNMKLMLLMLMKLMLMMJJMOOLLLNMMNEEEEEEEEEFFFNOMMNML	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPOOOOQOLNKMPPQNNLNOOOPQNPQPPNQPNORMRNOMRNQMPQQPNMNPNPOJJJJJJJKJJJJOQOQPQN	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:1:13308:3220:8618	99	18	29769796	60	27M1I47M	=	29769971	250	AGCCAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAAGCCAAAAAACCAACTCCGCAGGGCGGTGTTAAA	=>=@@BBC1C@5=A?@D>@:ACBDBE@@B44BBBBBBBBBBBBED.A-B-BD-B-/50D+,-42AE+A2D;2*==	MC:Z:75M	BD:Z:LLNKOONLNLNLLMLLLIJLNNKKKMLLNDEDDDEDEDDEEEMNKMNEEEEMMMNMNMMJMNNMJMNLOLLJLME	MD:Z:46C15A2T8	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPNQPNNQLQMNLQMQLOPQPNMNPNPOJJJJJJJJJJJJJOPNQPJJJJQOQPQRPMPPRPOMOQNQNNPOMJ	NM:i:4	MQ:i:60	AS:i:52	XS:i:26\n"
		"NS500817:604:HLY3VBGXC:3:23510:25620:14912	99	18	29769796	60	47M28S	=	29769965	244	AGCCAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAGCCAACCAACCCAAACCGGAGGGGGTGGTTAAAA	=>=@@BBCAC>@>A?@D>@E2C@DBAD@@BBBBBBBBAB..5-D..--;-A--B---D+2/4222214312*<==	MC:Z:75M	BD:Z:LLNKOONLNLNLLMLLLIJLNNKKKMLLNDEDDDEDEDDEEMNKMNMMMNMMIMNEMMJKMMMJJKOLNOLLMEE	MD:Z:47	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPNQPNNQLQMNLQMQLOPQPNMNPNPOJJJJJJJJJJJJOPNQPQOQPQONQPJQOPNPONNMMQMLRQOMJJ	NM:i:0	MQ:i:60	AS:i:47	XS:i:25\n"
		"NS500817:604:HLY3VBGXC:2:23204:14179:8481	163	18	29769799	60	75M	=	29770079	355	CAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATA	<=<=?>A>>=@>?D>@D@CACADC@@BBBBAABBBBBB9EC=B@ABBDA@BD5AD=EBE/E,AE<@?A>A@><;:	MC:Z:75M	BD:Z:LLNMMNLMKKLKLLIILNNKKKMLLMDDDDDDDDDDDDLNKMNMLNEMMMNMNMMJMNNLMNMNKLKKMOGEELK	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPNMQLQMNLQMPMOPPONLMONONIJJJJJJJJJJJOPNQPQNPJQOQPQRPMPPRPPOPQQMNNPONKJJRN	NM:i:0	MQ:i:60	AS:i:75	XS:i:27\n"
		"NS500817:604:HLY3VBGXC:4:22404:6543:13509	99	18	29769803	60	75M	=	29769986	258	CACTTAATTGTTGACACAGGAAAAAAAAAAAAAAGCCAACAACCCAACACCGCAGCGCTTTGTTAAAAATACTGG	=/<==<@>@D>@D3CACADDA@BBBBBAB4BBB.5E=B..@B?DDA-D.DD+,-4,+,5114?A?BBB85=.66=	MC:Z:75M	BD:Z:LLLNNMNLMMJJMOOKKKNLLMDDDDDDEDEDDDMMKLMMLNMMIMNMLLMJMNNNMMMLEMJJLNFFFMMLNMM	MD:Z:42A5T6A3G15	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMRNOMQMPLOPPOMMMOMPOJJJJJJJJJJJJOPNQPQNPQONQPQNMOPPRPPQPPNJQMQOMJJJSOQRQL	NM:i:4	MQ:i:60	AS:i:55	XS:i:26\n"
		"NS500817:604:HLY3VBGXC:4:12602:3682:13465	99	18	29769811	60	12M1I55M4S	=	29769998	262	TGTTGACACAGGAAAAAAAAAAAAAAAGCCAACAAACCAACTCCGCAGCGCTGTGTTAAAAATACTGGCGTT	==8<@>B@CADC@AAAAAABBABBBB.AED@..B.4D.;@.6E.,E-4-+,541@?A?BB770?-5E2-+==	MC:Z:75M	BD:Z:LLJJOPPLLLNMMNEDDDEDDDDDDDDLNJMMLKNDMLLNMNMMJMNNNMMMMJKJJLMEEELKLONNLMOJ	MD:Z:47A19	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMPQQPMLMOMONIIJIIIJJJJJJJOPNQPQNPJQOQPQRPMPPRPPQPPQMMMPOMJJJROQRQLOQQP	NM:i:2	MQ:i:60	AS:i:55	XS:i:25\n"
		"NS500817:604:HLY3VBGXC:4:12610:26254:5707	1123	18	29769811	55	26M42S	=	29769998	262	TGTTGACACAGGAAAAAAAAAAAAAACGCCCACAACCCAACCCCCCCGCCCGGGGTTAAAAAAAAGGG	=68<@>B@1ADC@AAAA6ABBABBB4.,E...D...DD........-@---@22@11?BBB778-4DD	MC:Z:75M	BD:Z:LLJJOPPLLLNMMNEDDDEDDDDDDDLMMJILKKNLMHLNMMIIIIIJMKIJKJJNJLMEEEEEEMMJ	MD:Z:26	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMPQQPMLMOMONIIJIIIJJJJJJQNPNNQMNPQONQPQONNNNNPPNNPNMMQPOMJJJJKJONM	NM:i:0	MQ:i:60	AS:i:26	XS:i:22\n"
		"NS500817:604:HLY3VBGXC:4:13408:10436:5182	163	18	29769825	60	75M	=	29769973	223	AGAAAAAAAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGG	<>8=>>?>?@@@CDCAACAAACCAADAED=D@EAEEAE>E>@0BABB;=DAEDADAAAEB@EDDAAEA?@C@;@?	MC:Z:75M	BD:Z:LLLNFEEDDDDDLMJLMLKMDLLLMLMLLILMMKLMLLIKJJLMEEELKLNMMMONLELNMMMJNMNPPOOMMMM	MD:Z:1A73	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOOOHIIIIIIINONPOPMPIPNQOPQPMPPRPOOPPQMMMPOMJJJRNQRQLPPRNJQRRQNNQORQPORLPON	NM:i:1	MQ:i:60	AS:i:73	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:1:12307:1747:13022	163	18	29769832	60	75M	=	29769992	235	AAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATT	<=;=>ABA?@B@@@CCAAC@DC?D@EAEEAD>E?A?BBAB?>DADDADAAAEB@EDDAAEA@AEDAED>>>=9;<	MC:Z:75M	BD:Z:LLEEFMNJLMLKMDLLLMLMLLILMMKLMLLIJIIKLDDELKLNMMMONLELNMMMIMLMOOMNNNNNOONMLLM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOJJHNOMPOPMOIQNPOPROLOPQONOPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNRQPOQLQPNQPMNOLN	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:3:21610:25350:5908	163	18	29769832	60	75M	=	29769955	198	AAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATT	<=;=>ABA?@B@@ACCAAC@DC?DAEAEEAD>E?A?BBA?>>DADDADAAAEB@EDD>AEA=AEDAED>>>=9;<	MC:Z:75M	BD:Z:LLEEFMNJLMLKMDLLLMLMLLILMMKLMLLIJIIKLDDELKLNMMMONLELNMMMIMLMOOMNNNNNOONMLLM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOJJHNOMPOPMOIQNPOPROLOPQONOPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNRQPOQLQPNQPMNOLN	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:4:23410:25485:10874	163	18	29769832	60	75M	=	29770036	279	AAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATT	<=;=>ABA?@B@@ACCAAC@DC?DAEAEEAD>E?A?BBAB?>DADDADAAAEB@EADAAEA@AE@AED>>>=9;<	MC:Z:75M	BD:Z:LLEEFMNJLMLKMDLLLMLMLLILMMKLMLLIJIIKLDDELKLNMMMONLELNMMMIMLMOOMNNNNNOONMLLM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOJJHNOMPOPMOIQNPOPROLOPQONOPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNRQPOQLQPNQPMNOLN	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:1:12201:18063:20373	99	18	29769839	60	75M	=	29770002	238	CAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCAT	==<?@?@BCAAC@DC?DADAEDAE?E?@?BBBB@?DAEDADAAAEB@EDDAAEA@AEDAED?>@A?@A?C@A?=;	MC:Z:75M	BD:Z:LLNMNOFMMMNMNMMILMNKLMLLIJIILLEDDKKKNLLMONLELNMMMIMLMOOMMMMMMNMLMMMNMMNNNNM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPQNPJPNPOPQOLOPQONOPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOMNOLNOROPPRR	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:3:23606:3979:3985	163	18	29769839	60	75M	=	29769962	198	CAACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCAT	<=;=>>?AA@@B?DC?DAD@DD@D>E?0?BAAB@?DAEC?C@AADB@EDDAAEA@AEDAED?>@A?@A>C@A?=;	MC:Z:75M	BD:Z:LLNMMNELLLMLMLLILMMKLMLLIJIIKLDDDKJKMLLMONLELNMMMIMLMOOMMMMMMNMMNMMNMNNNNNM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPQLOIPNPOPQOMOOQOONOOQLLLPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPPNQPLNPMNOROPPRR	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:2:22209:7220:12764	1123	18	29769840	60	75M	=	29770032	265	AACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATC	==<?@?BBAAC@D1=AAA@AE@7?E?A>BBBB@?DAEDADAAAEB@EDDAAEA@AEDAED?>@A?@A?DA@A<;?	MC:Z:73M	BD:Z:LLMLPFNMMNMNMMJLMMLLMLLIJIIKMDEDKJLMMLLONLELNMMMIMLMOOMMMMMMNMLMLMNMMMPNNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQNPJQNPOPQOLOORONNPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLOOLNOQOQPRRQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:4:11507:25875:3744	99	18	29769840	60	75M	=	29770030	265	AACAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATC	==<?@?BBAAC@DC?DAD@EE@E?E?A>BBBB@?DAEDADAAAEB@EDDAAEA@AEDAED?>@A?@A?DACA<;@	MC:Z:75M	BD:Z:LLMLPFNMMNMNMMJLMMLLMLLIJIIKMDEDKJLMMLLONLELNMMMIMLMOOMMMMMMNMLMLMNMMMPNNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQNPJQNPOPQOLOORONNPPQMMMPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLOOLNOQOQPRRQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:2:21305:8003:6690	99	18	29769841	60	75M	=	29769976	210	ACAAACCAACTCCGCAGAGCTGTATTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCA	==/?@>?@@@@4C/AA;@DEA@4,@A=ABBB@?DAE3A.AAAEB@1BDAAEA@A@DA421>@A?@A?DBDC>:==	MC:Z:75M	BD:Z:LLLNGNNMNMNMMJMMMKMMLLILKLKLEDEKJKNLMLNNLELNMMMIMLMOOMMMMMMNMLMLLNMMMOPNMMN	MD:Z:23G51	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OONPJQOPOPQOLOOQPNNOPQMOLNOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNPLNOQNQQRRQQ	NM:i:1	MQ:i:60	AS:i:70	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:1:11311:1890:2971	163	18	29769842	60	75M	=	29769999	232	CAAACCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAA	<=;=>@>?A?CB>DAD@DD@D>D>@?BBBB?>DAEDAD@A@DB@DDD<AEA@AEDAED?>@A?@A?@BDD@<<==	MC:Z:75M	BD:Z:LLNENMMMLMLLILMMKLMLLIJIIKLDDDKJKMLLLNMLELNMMMIMLMOOMMMMMMNMLMLMNMMMOPPMMNN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPJONPOPQOLOORONNOPPLLMONLJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQPLNPLNPRNPQSRQQP	NM:i:0	MQ:i:60	AS:i:75	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:3:22506:16796:9689	99	18	29769845	60	75M	=	29769973	203	AACAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCT	==<?@A?CC?DAD@DD@D>E?@?BBBB>?DAEDADAAAEB@EDDAAEA@AEDAED>?@A?@A?DBEEB@D@>?@<	MC:Z:75M	BD:Z:LLMLPNOMMJMNNLMMLLJJIIKLDDDKKKNLLLOMLDKNMMMIMLMOOMMMMMMNMLMLLMLLLOOONNPNMNM	MD:Z:1C73	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQNPQROLOOQONNOPPLLMPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNORNPPRRRRPOPP	NM:i:1	MQ:i:60	AS:i:73	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:1:11202:22485:14670	1123	18	29769846	60	75M	=	29770094	323	CCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTG	=><?@>CB?DAD@DD@D>D?A>BBBB@>DAEDADAAAEB@EDDAAEA@AEDAED>?@A?@A?DBEEB@EA@A?<@	MC:Z:75M	BD:Z:LLMNOONMJMNNLMNLLIKIIKLDDDKJLMMLLNNKEKMMMMIMLMOOMMMMMMNMLMLLMLLLNOONNOPMNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQPQRPLOOQONNOOQLLLPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQOPPRRQRQOPPQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:1:11312:5110:7144	1123	18	29769846	60	75M	=	29770094	323	CCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTG	=><?@>CB?DAD@DD@D>D?=>BBBB@>DAE@0DAAAEB@EDBAAEA@AEDAED>?@A?@A?DBEEB@EA?A?<@	MC:Z:75M	BD:Z:LLMNOONMJMNNLMNLLIKIIKLDDDKJLMMLLNNKEKMMMMIMLMOOMMMMMMNMLMLLMLLLNOONNOPMNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQPQRPLOOQONNOOQLLLPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQOPPRRQRQOPPQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:12401:14102:13694	99	18	29769846	60	75M	=	29770094	323	CCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTG	=><?@>CB?DAD@DD@D>D?A>BBBB@>DAEDADAAAEB@EDDAAEA@AEDAED>?@A?@A?DBEEB@EA@A?<@	MC:Z:75M	BD:Z:LLMNOONMJMNNLMNLLIKIIKLDDDKJLMMLLNNKEKMMMMIMLMOOMMMMMMNMLMLLMLLLNOONNOPMNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQPQRPLOOQONNOOQLLLPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQOPPRRQRQOPPQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:12612:24627:14661	99	18	29769846	60	75M	=	29769988	217	CCAACTCCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTG	=><?@>CB?DAD@DD@D>D?A>BBBB@>DAEDADAAAEB@EDDAAEA@AEDAED>?@A?@A?DBEEB@EA@A?<@	MC:Z:75M	BD:Z:LLMNOONMJMNNLMNLLIKIIKLDDDKJLMMLLNNKEKMMMMIMLMOOMMMMMMNMLMLLMLLLNOONNOPMNMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQPQRPLOOQONNOOQLLLPOMJJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQOPPRRQRQOPPQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:23502:16064:18761	99	18	29769846	60	75M	=	29770045	274	CCAACTCCGCAGAGCTGTGTTGAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTG	=><?@>CB?DAD@DD@D>D?A7ABBB@>DAEDADAAAEB@EDDAAEA@AEDAED>?@A?@A?ABEEB@EA@A?<@	MC:Z:75M	BD:Z:LLMNOONMJMNNLMNLLIKIILNMDDKJLMMLLNNKEKMMMMIMLMOOMMMMMMNMLMLLMLLLNOONNOPMNMM	MD:Z:21A53	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQPQRPLOOQONNOOQLLLPQQOJJRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQOPPRRQRQOPPQ	NM:i:1	MQ:i:60	AS:i:70	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:22101:5413:10879	99	18	29769852	60	75M	=	29770029	252	CCGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGC	=>9@@B?CD@D>D>@>AAAB@>DAEDACAAAEB@EDDAAEA@AEDAED>?@A?@A?DBEEB@EBBEEAE@B@;@@	MC:Z:75M	BD:Z:LLJMPOMMNMMJKJJKLDEDKJKMLLLNNKEKMLMLILKMOOMMMMMMNMLMLLMLLLNNNMMNNNONNPQMMMN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPPRPONOOPLLLONMIIIRNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQRPOPPQRQOQQS	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:2:12209:22204:13954	99	18	29769853	60	75M	=	29770020	242	CGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCA	=:=?A>CC@D>D>@=AAAA@?CAEDAD@AAEB@EDDAAEA@AEDAED>?@A?@>?DBEEB@EBBEEAE0C@=?>=	MC:Z:75M	BD:Z:LLMNPMNNMMJKJJLLDDEKJKMLLLNMLDLMLLMHMKLOOMMMMMMNMLMLLMLLLNNNMMNNMONNPPOMMNN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRPOOOOPLLLONLJIIQNQRQLPPRNJQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQQQOPPQQQPQQSR	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:4:21611:13599:18883	161	18	29769853	60	75M	=	29770216	438	CGCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCA	<:<=?=BB=C=C=@>AAAA?>C@DC@DAAADA@EDD>ADA?@ED@ED>?@A?@A?DBEEB@E@BEEAE@CB=?@=	MC:Z:75M	BD:Z:LLMNOLMMLLIJIIKLDDDKJKMLLLNMKDKMLLLHLKLOOMMMMMMNMLMLLMLLLNNNMMNONONNPQOMMNN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRNNNOOPLLLOOLIIIRMPQQKOORNJQQRQMNQNQQPNQLPONQOLNOLNOQNPPSRQRPOQQQQQPQQSR	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:23105:6484:3315	163	18	29769854	60	75M	=	29770003	224	GCAGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAA	<>;>>AB>B=C=?>AAAA?>C@DC@D5AAEA?ED->1E@<1=DADD1?@A)@A=DBEEB@EBB4?AE<CC?A<==	MC:Z:75M	BD:Z:LLNNMMNLLIJIIKLDDDKJKMLLLNMKDKMLLLHLKLNOMMMMMMNMLMLLMLLLNNNMMNNNONNPPOOMNNN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OORPMNOOPLLLONMIIIQNPQPLOOQNJQQRQMNQNQQPNQLPONQOLNOLNOQNPPRSQQQOPQRQPPRQSRP	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:13111:9009:11794	163	18	29769856	60	75M	=	29769988	207	AGAGCTGTGTTAAAAATACTGGACTATCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATA	<>;>?=8<82?=@AAA?>C@4@@@@)>@B=D-?5=@A0@EC@4D>)@;?@A>?BE,B<@B@EEA@@@-<?C><:9	MC:Z:75M	BD:Z:LLLMOMMIJIIKLDDDKJKMLLLNMLKLMLLLHLKLNNLMMMMMNMLMLLMLLLNNNMMNNMNNNPPNNOPNNLK	MD:Z:25T49	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOOONOPLLLONLIJIQMPRPKOPQNKQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQRPOQPQRQOQRTRPRN	NM:i:1	MQ:i:60	AS:i:70	XS:i:21\n"
		"NS500817:604:HLY3VBGXC:3:21511:11891:1700	1187	18	29769856	60	75M	=	29769990	207	AGAGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATA	<>;:?6B1?=?/@4AA?>@@@2@C@AA0B@>-DAA@A@@E2=ED>?@A=@A?-BEEB@@BBE>AA/1=5;C.<;:	MC:Z:73M	BD:Z:LLLMOMMIJIIKLDDDKJKMLLLNMKDKMLLLHLKLNNLMMMMMNMLMLLMLLLNNNMMNNMNNNPPNNOPNNLK	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOOONOPLLLONLIJIQMPRPKOPQMIQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQRPOQPQRQOQRTRPRN	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:2:13312:13319:17587	163	18	29769858	60	75M	=	29769956	173	AGCTGTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATG	<>1;6;B<>=@@@A?>C@DC@C@@@EA@EDC>AEA@AECADC?)?A;@A?/BEEB@EBBEEAE/-DAE3A@<9;@	MC:Z:75M	BD:Z:LLNMNJKIIKLDDDKJKMLLLNMKDKMLLLHLKLNNLLLMMMNMLMLLMLLLNNNMMNNMNMMPPNNNOPPLKLN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPPOLLLONLIIIRMPQPLOOQNIPPRQMNQNQQPNQLPONQOLNOLNOQNPPRRQQPPPPRQPPRQSSQRNLP	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:13102:13541:6857	99	18	29769862	60	75M	=	29769989	202	GTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATC	=9<;=<@@AA?>C@DC@C@AADB@EDD@AEA@AEDAED>?@A?@A?DBEEB@EBBEEAEADDAEEBB@>?C@;;@	MC:Z:75M	BD:Z:LLKJLMNEEELKLNMLLNNKDKMLLLHLLLONLLMLMLMMLMLLMLLLNNNMMNNMNMMOOMMMNOOMLMPMMOM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMMPOMIIIQMPQPKPOQMJQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQRSRPRNMQLPPQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:23403:8907:13256	99	18	29769863	60	75M	=	29769928	140	TGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCA	==8<=?@@A?>C@DC@C@@AEA@EDDA@EA@AEDAED?>@A0@A?DBEEB@EBBEEAE/DDAEEBB@?=DB=:@=	MC:Z:75M	BD:Z:LLJJNNFEELKLNMMLNMLDKMLLLHLKMNOLLLMLMMLLMLLMLLLNNNMMNNMNMMOOMMMNNOMLMOOMOMN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMPOMJIIQMPQPKOPQMIQQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQQTRPRNLQMPPQQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:4:23609:1678:4094	147	18	29769863	60	6S69M	=	29769864	-68	CCGATCTGTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATAT	?;?<>DA@8=AAAA@=?E@DE?E@BAE@ADCEB-D@ABDEADD>?B?@A)>D@DD?@D@@DD>C=BB?BA;:*:/	MC:Z:68M7S	BD:Z:KOMONPMNNMFFEMLMPNMMOPMENONOMJMMNNMLMLLLLLKJKLJKLLILLMMNNLKLMMMNMLMMNMNLKLL	MD:Z:69	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:NNQPPQNQNPJJKNLPPPQMQPOJOQPPPMNOPQQRPQMPQOQNRMNRLNLQOROPPPMOPOPPPMOQRPLLNOO	NM:i:0	MQ:i:60	AS:i:69	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:4:23609:1678:4094	99	18	29769864	60	68M7S	=	29769863	68	GTTAAAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATAGATCGG	=99<@>>@?>C@DC@C@@@E7?EADAADA@AED8ED>?@A/@A?DBEE@@E-BEEAE:@DAEE@-@>@?D?<<:?	MC:Z:6S69M	BD:Z:LLJLOFFELKLNMMMNMKEKMLLLHLKLONMLLLMLNLKMLLMLLLNNNMMNNMNMMOOMMMNNNMLMLONOMOK	MD:Z:68	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPOMJJIQMPQPKOORMIPQRQMNQNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQQSSPRNLNQPPQNN	NM:i:0	MQ:i:60	AS:i:68	XS:i:20\n"
		"NS500817:604:HLY3VBGXC:3:22605:25184:12691	99	18	29769868	60	74M	=	29769988	195	AAAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGA	==<?=<>>@C@C@@@DA4@DD@A@A@A@@A93?>=A?@A?@BEE.@EBBEB5EADAA@EBB4?@E2A@EA>2<<	MC:Z:75M	BD:Z:LLEENLMNMMMONLEKMLMLHLKLNNLLMLMLMLLLLKLLLLNNNMMNNMNMMOOMMMNNNLKLNNNPNOMMMO	MD:Z:74	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOJJRNQQPKOOQMIPQQPLNQNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQQSRPRNMPLPPQRRNQQ	NM:i:0	MQ:i:60	AS:i:74	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:23406:14666:7551	99	18	29769869	60	74M	=	29770057	263	AAATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAG	==6<=A?CC@@@@@DA?DCD>@EA@AECAED>?@A?@A?DBEEA@EB@EEAE9DDAEEBB@?@EDA@EB>?>;@	MC:Z:75M	BD:Z:LLELMMOMMMONLELMLLMHLKLNNLLLMLMMLKMKLLKLLNNNMMNNMNMMOOMMMNNNLKLNMNPNONMMOM	MD:Z:74	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOJRNQRPKOOQMIPPRPLMQNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQQSRPRNLQLPPQQSNQQO	NM:i:0	MQ:i:60	AS:i:74	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:11609:3802:11536	99	18	29769870	60	75M	=	29769985	190	AATACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGAT	==9<@>CB@C@@@DA?DCCAADA@AED@ED>?@A?@A?DBEEB@EBBEEAEADDAEEBB@?@EDA@EB@@C=?<;	MC:Z:75M	BD:Z:LLLKNONMMONLELNLLLILKLNNLLLLMLNLKLLKMKKLNNNMMNNMNMMOOMMMNNNLKLNMMPNONNOOMLO	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OORNQRQKOOQMIPPQQLMPNQQPNQLPONQOLNOLNOQNPPRRQQPOPPQQPOQQSRPRNLPMPPQQRORQOOP	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:22610:3899:12350	163	18	29769872	60	75M	=	29769969	172	TACTGGACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCA	<:;;>@>A>>?C@?DCC@@D@?@D@AED>?>@?@A?DBDEA?EBAE?AE/DDAEEBB@?@EDA@EB@AD@C=:@=	MC:Z:75M	BD:Z:LLLNNMMNMKDKMLLLHLKLNNLLLLLLMLKLKKLKKKMNNMMNNMNMMOOMMMNNNLKLNMMPNONNNQOLOMN	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQROKOOQMIPPQQLMPMQPOMQKONNQOLNOLNOQNPPRRQQPOPPQQPOQQSRPRNMPLQPQRSNQRPOPQQ	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:1:12203:14253:3380	163	18	29769877	60	75M	=	29770177	375	GACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGC	<=;;;=?>=@B@?@D@?@DC@DC>>@A?@A>CBEEB@EA?DD>E>DDAEEBB@?@EDA@EB@AEAEA=DA@>>?@	MC:Z:75M	BD:Z:LLONMELMLLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMNMNMMOOMMMNNNLKLNMMOMNMMNPNMPNPPEMML	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRLIPPQPLMPMQPOMPLONMQNKMOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPRQROQQPPPQRQJONO	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:2:21206:22463:9954	1187	18	29769877	60	75M	=	29770157	355	GACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGC	<=;;;=B?=CBB?2@@?@@C@D52>@>?@A>CBEEB@@ABDDAE/DD<@E@B@?@E2A@EB0A4AEA@DA@>=?@	MC:Z:75M	BD:Z:LLONMELMLLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMNMNMMOOMMMNNNLKLNMMOMNMMNPNMPNPPEMML	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRLIPPQPLMPMQPOMPLONMQNKMOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPRQROQQPPPQRQJONO	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:12605:2634:6064	163	18	29769877	60	75M	=	29770157	355	GACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGC	<=;;;=B?=CBB?@D@?@DC@DC>>@A?@A>CBEEB@EABDDAE>DDAEEBB@=@EDA@EB@AEAEA@DA@>>?@	MC:Z:75M	BD:Z:LLONMELMLLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMNMNMMOOMMMNNNLKLNMMOMNMMNPNMPNPPEMML	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRLIPPQPLMPMQPOMPLONMQNKMOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPRQROQQPPPQRQJONO	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:3:21503:12374:7349	1187	18	29769877	60	73M	=	29770157	355	GACTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAG	<=;;;=1-=CBB?2D@11DC@DC>>@A?@=>CBEE@@EABDDAE@DDAEEBB@=@4DA@EB@AE@AA@?A@>?	MC:Z:75M	BD:Z:LLONMELMLLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMNMNMMOOMMMNNNLKLNMMOMNMMNPNMPNNNEM	MD:Z:73	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOPRLIPPQPLMPMQPOMPLONMQNKMOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPRQROQQPPPQQPJO	NM:i:0	MQ:i:60	AS:i:73	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:4:12610:1610:19447	1187	18	29769877	60	75M	=	29770157	355	GCCTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGC	<-<7;=1>1CBB?@D0=@DC@D2>=@A?@A>@BEA-@E--4DA@/D-A@EBB@?@ED/@EB@AE/4/=?A@>5?@	MC:Z:75M	BD:Z:LLKMMELMLLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMNMNMMOOMMMNNNLKLNMMOMNMMNPNMPNPPEMML	MD:Z:1A73	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OONQLIPPQPLMPMQPOMPLONMQNKMOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPRQROQQPPPQRQJONO	NM:i:1	MQ:i:60	AS:i:73	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:1:12201:18355:2938	163	18	29769879	60	75M	=	29770056	252	CTTTCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGCTA	<;8:>>=BAB>?C@1@DC@DC>>?@?@A?DADEB@EBBDE@D/DCAEEBB@?@EDA@EB@A;AE@@EBAAC@?<:	MC:Z:75M	BD:Z:LLLEMNMLLHLKLNNLLLLLLMLKLKKLKKKMMMLLMMLNMMOOMMMNNNLKLNMMOMNMMMONMPNOOGOMLMM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OONJOPQPLMPMPPPMPKOOMPNLMNKNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPQQSNQROOQRQPKPNOPO	NM:i:0	MQ:i:60	AS:i:75	XS:i:19\n"
		"NS500817:604:HLY3VBGXC:4:13510:18834:6681	99	18	29769882	60	75M	=	29770090	283	TCATCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGCTATAA	==<<@AB?@D@?@DC@DC>?@@?@A?DAEEB@EBBEEAEADDAEEBB@?@EDA@EB@AEAEA@EBBBEDD?;::=	MC:Z:75M	BD:Z:LLNMONJMLMOOMMMLLLNLKLKKLKKKNMNLLMNLNLLOOMMMNNNLKLNMMOMNMMMOMLOMNOFNNMOMLKM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOQRQMNPMPPOMPKOOMPNLNOLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPQQRNQQOOPRQPJONPQOLNM	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:22301:22834:17268	163	18	29769885	60	75M	=	29770003	193	TCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGCTATAAATT	<=<><=B>=?CB?DC>>?@>?@>CAEEB@EAAEEAEADCADDBB>?@EDA@EB@AEAEA@EBBBEDEA>?=><;<	MC:Z:75M	BD:Z:LLMINLMNNLLLLLLMLKLKKLKKKMMMLLMMLMLLNNLMMNNNLKLNMMOMNMMMOMLOMNNFNNMNNNMMELM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMNOMPPOMPKONNPNKMOKMNQMOORRQQPOPPQQPOQQSRPRNLPLPPQQRNQQOOQQQQJOOPPOMOMJRN	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:23201:4271:15898	1187	18	29769885	60	75M	=	29770075	265	TCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGCTATAAATT	<=<>7=B><?CB?DC>=?@>?@>CAEEB@EA4EEAE/@CADDB->?@EDA@EB@AEAE/@EBB-EDE>)?=><;<	MC:Z:75M	BD:Z:LLMINLMNNLLLLLLMLKLKKLKKKMMMLLMMLMLLNNLMMNNNLKLNMMOMNMMMOMLOMNNFNNMNNNMMELM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMNOMPPOMPKONNPNKMOKMNQMOORRQQPOPPQQPOQQSRPRNLPLPPQQRNQQOOQQQQJOOPPOMOMJRN	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n"
		"NS500817:604:HLY3VBGXC:1:23312:19902:7696	99	18	29769885	60	75M	=	29769980	170	TCCCTTGATTGGAGGTATTATTACAGCATCAAGCTGACCTGCAATATGGATCATTGAGATCAAAGGCTATAAATT	===@>>C>?@DC@DC>>?@?@@?DBEEA@EBBEEAEADDAEEBB@?@EDA@EB@AEAEA@EBBBEDEA>?=><;<	MC:Z:75M	BD:Z:LLMIOMNOOMMMMMMMLKMKKLKKKMMMMLNMLMMLONLMMNNNLKLNMMOMNMMMOMLOMNNEMNMNNMMMELM	MD:Z:75	PG:Z:MarkDuplicates	RG:Z:10565081	BI:Z:OOMNQNQPOMPKONMPOKMNLNOQNPPRRQQPOPPQQPOQQSRPRNLPLPPQQRNQQOOPQQPKONOPOMOMJRN	NM:i:0	MQ:i:60	AS:i:75	XS:i:0\n";
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace std;
    using namespace varoom;

    istringstream in(sam_text);
    sam_reader S(in);
    BOOST_CHECK_EQUAL(S.more(), true);
    size_t n = 0;
    for(; S.more(); ++S)
    {
        ++n;
    }
    BOOST_CHECK_EQUAL(n, 48);
}

BOOST_AUTO_TEST_CASE( test2 )
{
    using namespace std;
    using namespace varoom;

    istringstream in(sam_text);
    sam_reader S(in);
    BOOST_CHECK_EQUAL(S.more(), true);
    size_t n = 0;
    for(; S.more(); ++S)
    {
        ++n;
        switch (n)
        {
            case 2:
            {
                const sam_alignment& aln = *S;
                BOOST_CHECK_EQUAL(aln.name, "NS500817:604:HLY3VBGXC:1:13308:3220:8618");
                BOOST_CHECK_EQUAL(aln.flags, 99);
                BOOST_CHECK_EQUAL(aln.chr, "18");
                BOOST_CHECK_EQUAL(aln.pos, 29769796);
                BOOST_CHECK_EQUAL(aln.mapq, 60);
                BOOST_CHECK_EQUAL(aln.cigar, "27M1I47M");
                BOOST_CHECK_EQUAL(aln.mate_chr, "=");
                BOOST_CHECK_EQUAL(aln.mate_pos, 29769971);
                BOOST_CHECK_EQUAL(aln.tlen, 250);
                BOOST_CHECK_EQUAL(aln.seq, "AGCCAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAAGCCAAAAAACCAACTCCGCAGGGCGGTGTTAAA");
                BOOST_CHECK_EQUAL(aln.qual, "=>=@@BBC1C@5=A?@D>@:ACBDBE@@B44BBBBBBBBBBBBED.A-B-BD-B-/50D+,-42AE+A2D;2*==");
                BOOST_CHECK_EQUAL(sam_flags::is_paired(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_proper_pair(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_unmapped(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_mate_unmapped(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_reverse(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_mate_reverse(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_read_1(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_read_2(aln.flags), false);
                break;
            }
            case 4:
            {
                const sam_alignment& aln = *S;
                BOOST_CHECK_EQUAL(aln.name, "NS500817:604:HLY3VBGXC:2:23204:14179:8481");
                BOOST_CHECK_EQUAL(aln.flags, 163);
                BOOST_CHECK_EQUAL(aln.chr, "18");
                BOOST_CHECK_EQUAL(aln.pos, 29769799);
                BOOST_CHECK_EQUAL(aln.mapq, 60);
                BOOST_CHECK_EQUAL(aln.cigar, "75M");
                BOOST_CHECK_EQUAL(aln.mate_chr, "=");
                BOOST_CHECK_EQUAL(aln.mate_pos, 29770079);
                BOOST_CHECK_EQUAL(aln.tlen, 355);
                BOOST_CHECK_EQUAL(aln.seq, "CAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAGCCAACAAACCAACTCCGCAGAGCTGTGTTAAAAATA");
                BOOST_CHECK_EQUAL(aln.qual, "<=<=?>A>>=@>?D>@D@CACADC@@BBBBAABBBBBB9EC=B@ABBDA@BD5AD=EBE/E,AE<@?A>A@><;:");
                BOOST_CHECK_EQUAL(sam_flags::is_paired(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_proper_pair(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_unmapped(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_mate_unmapped(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_reverse(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_mate_reverse(aln.flags), true);
                BOOST_CHECK_EQUAL(sam_flags::is_read_1(aln.flags), false);
                BOOST_CHECK_EQUAL(sam_flags::is_read_2(aln.flags), true);
                break;
            }
        }
    }
    BOOST_CHECK_EQUAL(n, 48);
}
