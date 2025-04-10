#!/bin/bash

rm Production_MD/*

perl write_mdp.pl md.mdp; mv md_* ./Production_MD


