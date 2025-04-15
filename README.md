# Stark Interference Simulation for BaF Molecules

This repository contains simulation code and logs related to studying systematic error signals (W) in Stark interference measurements using BaF molecules. The project investigates why a linear relationship arises between the extracted W signal and the product of laser detuning and non-reversing (NR) electric field strength.

## Motivation

Experiments using <sup>138</sup>BaF molecules showed unexpected linear behavior between systematic error signals and experimental imperfections. This simulation aims to reproduce and understand that effect in a simplified three-level quantum system.

## Model Overview

The system consists of:
- Two nearly degenerate ground states (odd/even parity)
- One excited state (even parity)
- Time-dependent electric fields (reversing, depletion laser, and non-reversing)
- Population transfer computed as molecules pass through the interaction region

## Goals

- Reproduce the experimental behavior of W as a function of laser detuning × NR field
- Understand underlying physics causing the systematic linear dependence

## Contents

- `simulation_log.html`: Log of simulation runs and parameter sweeps
- `code/`: (If included) Python or MATLAB scripts used for numerical simulation
- `README.md`: Project overview and instructions

