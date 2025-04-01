

# Multi-objective TSP with Profits and Passengers (MO-TSPPP)  
**Exact, Heuristic, and Evolutionary Approaches**  

This repository contains the implementation, benchmark instances, and experimental data for the paper:  

> **"Multi-objective Traveling Salesman Problem with Profits and Passengers: Exact, Heuristic, and Evolutionary Approaches"**  
> *Presented at the Genetic and Evolutionary Computation Conference (GECCO 2025)*  

## 📄 Abstract  

This paper addresses the Multi-objective Traveling Salesman Problem with Profits and Passengers (MoTSPPP), an extension of the Traveling Salesman Problem with Profits, to model ridesharing systems in multi-objective contexts. This problem, which considers three objective functions (minimizing travel cost, minimizing travel time, and maximizing the bonuses collected), has not been studied before. We implement a mathematical formulation to compute Pareto-optimal solutions, two non-evolutionary heuristics, and two evolutionary algorithms (NSGA-II and MOEA/D). Non-evolutionary heuristics solve the salesman route, bonus collection, and passenger assignments sequentially, whereas NSGA-II and MOEA/D combine these three decision levels via genetic operators. Experiments on 84 instances with up to 200 vertices assess the difficulty of optimally solving the problem and the performance of the algorithms.

## 📂 Repository Structure  
```
.
├── /Code/               # Source code (C/C++)
├── /Instances/          # Benchmark datasets
├── /results/            # Raw experimental and processed results
```

## 👥 Authors  
| Name | Affiliation | Contact |  
|------|-------------|---------|  
| **Juvenal B. A. Silva** | Federal University of Bahia, Institute of Computing | [juvenal.bruno@ufba.br](mailto:juvenal.bruno@ufba.br) |  
| **Tarcisio A. P. Brito** | Federal University of Bahia, Institute of Computing | [tarcisio.brito@ufba.br](mailto:tarcisio.brito@ufba.br) |  
| **Ricardo A. Rios** | Federal University of Bahia, Institute of Computing | [ricardoar@ufba.br](mailto:ricardoar@ufba.br) |  
| **Tatiane N. Rios** | Federal University of Bahia, Institute of Computing | [tatiane.nogueira@ufba.br](mailto:tatiane.nogueira@ufba.br) |  
| **Islame F. C. Fernandes** | Federal University of Bahia, Institute of Computing | [islame.felipe@ufba.br](mailto:islame.felipe@ufba.br) |  


## 📜 License  
This project is licensed under the **ACM License**.  

## 📧 Contact  
For questions, contact the lead author: [juvenal.bruno@ufba.br](mailto:juvenal.bruno@ufba.br).  

---

