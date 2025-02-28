#Rscript scripts/SimulateTrees.R 5 1 1 0.9 100 10 > out1

#export JAVA_HOME=`/usr/libexec/java_home -v 17.0.13`
#
#for i in {1..10}
#do
#  Rscript scripts/MakeSimulateDNAxml.R simulations/5_1_1_0.9_100_trees/$i xml_templates/simulateDNA.xml
#  java -jar beastFX.jar simulations/5_1_1_0.9_100_trees/$i/simulateDNA$i'.xml' > simulations/5_1_1_0.9_100_trees/$i/sim_seq_out
#    Rscript scripts/MakeRunxml.R simulations/5_1_1_0.9_100_trees/$i xml_templates/estimateTree.xml 0.5
#done


#Rscript scripts/SimulateTrees.R 5 0.9 1 0.9 100 100 > out2


#export JAVA_HOME=`/usr/libexec/java_home -v 17.0.13`
#
#for i in {1..100}
#do
#  Rscript scripts/MakeSimulateDNAxml.R simulations/5_0.9_1_0.9_100_trees/$i xml_templates/simulateDNA.xml
#  java -jar beastFX.jar simulations/5_0.9_1_0.9_100_trees/$i/simulateDNA$i'.xml' > simulations/5_0.9_1_0.9_100_trees/$i/sim_seq_out
#    Rscript scripts/MakeRunxml.R simulations/5_0.9_1_0.9_100_trees/$i xml_templates/estimateTree.xml 0.5
#done

#Rscript scripts/SimulateTrees.R 0.5 0.1 0.1 0.1 100 100

#sim_par="0.5_0.1_0.1_0.1_100_trees"

#for i in {1..100}
#do
#  Rscript scripts/MakeSimulateDNAxml.R simulations/$sim_par/$i xml_templates/simulateDNA.xml
#  java -jar beastFX.jar simulations/$sim_par/$i/simulateDNA$i'.xml' > simulations/$sim_par/$i/sim_seq_out
#    Rscript scripts/MakeRunxml.R simulations/$sim_par/$i xml_templates/estimateTreeIIBD.xml 0.5 10
#done

#cd simulations/$sim_par/
#
#for i in {1..100}
#do
#    cd $i
#    java -jar ../../../stratigraphic-ranges.jar -version_file /Users/aga122/Git/stratigraphic-ranges/version.xml $i_prop0.5.xml > estimate_tree_out
#    cd ../
#done


#Rscript scripts/SimulateTrees.R 0.2 0.001 0.05 0.9 100 100 > out4
#
#sim_par="0.2_0.001_0.05_0.9_100_trees"
#
#for i in {1..100}
#do
#    Rscript scripts/MakeSimulateDNAxml.R simulations/$sim_par/$i xml_templates/simulateDNA.xml
#    java -jar beastFX.jar simulations/$sim_par/$i/simulateDNA$i'.xml' > simulations/$sim_par/$i/sim_seq_out
#    Rscript scripts/MakeRunxml.R simulations/$sim_par/$i xml_templates/estimateTreeIIBD.xml 0.5 20
#done


#cd simulations/$sim_par/

#for i in {1..100}
#do
#    cd $i
#    java -jar ../../../stratigraphic-ranges.jar -version_file /Users/aga122/Git/stratigraphic-ranges/version.xml $i'_prop0.5.xml' > estimate_tree_out
#    cd ../
#done

#for i in {1..100}
#do
#    Rscript ../../scripts/CompareEstimatedTransmissionDates.R $i 0.5
#done

#age simulations
#Rscript scripts/SimulateTrees.R 0.3 0.001 0.05 0.9 20 100 > out_age

sim_par="0.3_0.001_0.05_0.9_20_trees"
#
#for i in {1..100}
#do
#    Rscript scripts/MakeSimulateDNAxml.R simulations/$sim_par/$i xml_templates/simulateDNA.xml
#    java -jar beastFX.jar simulations/$sim_par/$i/simulateDNA$i'.xml' > simulations/$sim_par/$i/sim_seq_out
#    Rscript scripts/MakeRunxml.R simulations/$sim_par/$i xml_templates/estimateTreeIIBD.xml 0.5 20
#done

cd simulations/$sim_par/

for i in {1..100}
do
#    cd $i
#    java -jar ../../../stratigraphic-ranges.jar -version_file /Users/aga122/Git/stratigraphic-ranges/version.xml $i'_prop0.5.xml' > estimate_tree_out
#    cd ../
        Rscript ../../scripts/CompareEstimatedTransmissionDates.R $i 0.5
done

cd ../../

