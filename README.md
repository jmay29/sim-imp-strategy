# sim-imp-strategy
Code for: A real data-based simulation procedure to select an imputation strategy for mixed-type trait data.

## **Code for: A real data-based simulation procedure to select an imputation strategy for mixed-type trait data.**

### **Description**
This project entails a real data-based simulation study that compares the performance of imputation methods using a mixed-type trait dataset. Missing completely at random (MCAR), missing at random (MAR), and missing not at random (MNAR) scenarios were simulated. The accuracy of different imputation methods were evaluated with and without phylogeny in the form of gene trees.

### **Data**
Trait data for this project were obtained from:

> Meiri, S. (2018). Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecology and Biogeography, 27(10), 1168–1172. https://doi.org/10.1111/geb.12773

> Meiri, S. (2019). Data from: Traits of lizards of the world: Variation around a successful evolutionary design. Dryad Dataset. https://doi.org/10.5061/dryad.f6t39kj

COI sequence data for this project were obtained from the Barcode of Life Data System (BOLD):

> Ratnasingham, S., & Hebert, P. D. N. (2007). bold: The Barcode of Life Data System (http://www.barcodinglife.org). Molecular Ecology Notes, 7(3), 355–364. https://doi.org/10.1111/j.1471-8286.2007.01678.x. dx.doi.org/10.5883/DS-IMPMIX2

RAG1 and c-mos sequence data for this project were obtained from:

> Pyron, R. A., Burbrink, F. T., & Wiens, J. J. (2013). A phylogeny and revised classification of Squamata, including 4161 species of lizards and snakes. BMC Evolutionary Biology, 13(1), 93. https://doi.org/10.1186/1471-2148-13-93

> Pyron, R. A., Burbrink, F. T., & Wiens, J. J. (2013). Data from: A phylogeny and revised 	classification of Squamata, including 4161 species of lizards and snakes. Dryad Dataset. 	https://doi.org/10.5061/dryad.82h0m

### **Glossary**

* **Imputation**
  * MCAR = Missing completely at random
  * MAR = Missing at random
  * MNAR = Missing not at random
  * KNN = K-nearest neighbour
  * RF = Random forest
  * MICE = Multivariate imputation by chained equations
  * PMM = Predictive mean matching
  * LR = Logistic regression


* **Performance measures**
  * MSE = Mean squared error
  * PFC = Proportion falsely classified
  * ER = Error ratio (error rate without phylogeny/error rate with phylogeny)

* **Traits**
  * AT = Activity time
  * IE = Insular/endemic
  * LC = Largest clutch
  * SC = Smallest clutch
  * SVL = Snout-vent length
  
  
* **Phylogenetics**
  * BM = Brownian motion
  * K = Blomberg's K (Blomberg et al. 2003 Evolution)
  * D = Fritz and Purvis' D (2010 Conservation Biology)


* **Genes**
  * COI = Cytochrome c oxidase I gene
  * CMOS = Oocyte maturation factor Mos gene
  * RAG1 = Recombination activating gene 1
