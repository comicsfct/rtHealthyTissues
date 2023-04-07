# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 01:39:08 2023

@author: paulo
"""

tissue_colors = {
'adipose-subcutaneous':'#DCD8D8', # light gray
'adipose-visceral':'#BDBBBB', # dark gray
'artery-aorta':'#FF5555', #red
'artery-coronary':'#FF6666', #red
'artery-tibial':'#FF8888', #red
'brain-cerebellum':'#FF9911', 
'brain-cortex':'#FF9933' , #orange
'breast-mammary-tissue':'#EC211A', # dark orange
'colon-sigmoid':'#D2D200', # dark yellow
'esophagus-gastroesophageal-junction':'#FFE666', #yellow
'esophagus-mucosa': '#FFE666', #yellow
'esophagus-muscularis':'#FFE666', #yellow
'heart-left-ventricle':'#91D31A', #green
'heart-atrial-appendage':'#91D31A', #green
'lung':'#48BF91', #violet
'liver':'#3ABDB3', #green sea
'muscle-skeletal': '#4666B4',
'nerve-tibial':'#1682B4',
'pituitary':'#4682B4', #cyan
'skin-sun-exposed':'#2682B4',
'skin-not-sun-exposed':'#2682B4',
'testis':'#BF80FF', #pink
'thyroid':'#FF80EA', #blue
}

tissue_abbrev = {
'adipose-subcutaneous':['Adipose-Subcut', '#FFE14A'],
'adipose-visceral': ['Adipose-Visceral', '#E3C942'],
'artery-aorta':['Artery-Aorta', '#FF6666'],
'artery-coronary':['Artery-Coronary', '#FF7777'],
'artery-tibial': ['Artery-Tibial', '#FF8888'],
'brain-cerebellum':['Brain-Cerebellum', '#76A2E3'], 
'brain-cortex': ['Brain-Cortex', '#6E97D4'],
'breast-mammary-tissue': ['Breast-Mammary', '#E4B5ED'],
'colon-sigmoid':['Colon-Sigmoid', '#B56D37'],
'esophagus-gastroesophageal-junction': ['Esophagus-Gastro', '#FF9A4D'],
'esophagus-mucosa': ['Esophagus-Mucosa', '#E68B45'],
'esophagus-muscularis': ['Esophagus-Muscle', '#C7783C'],
'heart-atrial-appendage': ['Heart-Atrial-App', '#BF3A3A'],
'heart-left-ventricle': ['Heart-L-Ventricle', '#AB3434'],
'lung': ['Lung', '#74C296'],
'liver': ['Liver', '#75825B'],
'muscle-skeletal': ['Muscle-Skeletal', '#FF9159'],
'nerve-tibial': ['Nerve', '#6084BA'],
'pituitary': ['Pituitary', '#5575A3'],
'skin-sun-exposed': ['Skin-Not-Exposed', '#E1FF4A'],
'skin-not-sun-exposed':[ 'Skin-Exposed', '#D6FF4F'],
'testis': ['Testis', '#7DD43B'],
'thyroid': ['Thyroid', '#876FAB'],
}

tissue_acronyms = dict([[key, value[0]] for key,value in zip(tissue_abbrev.keys(), tissue_abbrev.values())])
tissue_colors = dict([[key, value[1]] for key,value in zip(tissue_abbrev.keys(), tissue_abbrev.values())]).values()