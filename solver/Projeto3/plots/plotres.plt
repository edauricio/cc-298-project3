#!/usr/bin/bash --persist

plot "~/Workspace/C++/CC-298/cc-298-project3/solver/Projeto3/plots/RES.data" u 1:2 w lines title 'Continuity', \
     "~/Workspace/C++/CC-298/cc-298-project3/solver/Projeto3/plots/RES.data" u 1:3 w lines title 'X-Momentum', \
     "~/Workspace/C++/CC-298/cc-298-project3/solver/Projeto3/plots/RES.data" u 1:4 w lines title 'Y-Momentum', \
     "~/Workspace/C++/CC-298/cc-298-project3/solver/Projeto3/plots/RES.data" u 1:5 w lines title 'Energy'
