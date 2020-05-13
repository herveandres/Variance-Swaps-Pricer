cmake_minimum_required(VERSION 3.9.1)

set (CMAKE_CXX_STANDARD 11)

# set the project name
project(VarianceSwapsPricer)

# add the executable
add_executable(VarianceSwapsPricer 
                main.cpp 
                Model.cpp Model.h
                PathSimulator.cpp PathSimulator.h
                HestonLogSpotPathSimulator.cpp HestonLogSpotPathSimulator.h
                HestonVariancePathSimulator.cpp HestonVariancePathSimulator.h
                VarianceSwap.cpp VarianceSwap.h
                VarianceSwapsPricer.cpp VarianceSwapsPricer.h
                VarianceSwapsHestonAnalyticalPricer.cpp VarianceSwapsHestonAnalyticalPricer.h 
                VarianceSwapsHestonMonteCarloPricer.cpp VarianceSwapsHestonMonteCarloPricer.h
                MathFunctions.cpp MathFunctions.h) 