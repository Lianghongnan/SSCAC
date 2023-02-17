% %%
clc,clear;
folder_now = pwd;
addpath([folder_now, '\funs']);
runtimes = 10;
load ('gal.mat', 'new_Gal');
DATA=new_Gal;
result=main(DATA, runtimes)




