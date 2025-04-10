#!/bin/bash

rm EM/*
rm NVT/*
rm NPT/*
rm Production_MD/*

perl write_mdp.pl em_steep.mdp; mv em_steep_* ./EM
perl write_mdp.pl nvt.mdp; mv nvt_* ./NVT
perl write_mdp.pl npt.mdp; mv npt_* ./NPT
perl write_mdp.pl md.mdp; mv md_* ./Production_MD


