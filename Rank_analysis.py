import os
import sys
import itertools
from copy import deepcopy
import random
import queue
from collections import Counter,defaultdict
import pandas as pd
import numpy as np

def get_combinations(target_list=None,r=2):
    """
    return a iterable containing the combinations
    """
    if target_list:
        index_list=target_list
    else:
        index_list=[i for i in range(23)]
    comb_list=itertools.combinations(index_list,r)
    return comb_list

def read_extract(file_name,work_dir,cal_control,shuffle=False):
    """
    input: file_name,work_dir,control or diseases,shuffle
    return: numpy.array
    """
    file_name1=os.path.join(work_dir,file_name)
    df=pd.read_csv(file_name1)
    #### accodring to the index
    #### control:0, case:1
    if cal_control:
        label=0
    else:
        label=1
    if shuffle:
        label_list=list(df.iloc[:,-1])
        random.shuffle(label_list)
        ### resign the label column
        df.iloc[:,-1]=label_list
    ### the last columns is the label columns
    index=df.iloc[:,-1]==label
    tumor_df=df.loc[index,]
    ### columns to remove(first:names,last:labels)
    column_indices=[0,len(tumor_df.columns)-1]
    drop_names=tumor_df.columns[column_indices]
    tumor_df_M=tumor_df.drop(drop_names,axis=1)
    ny_df=tumor_df_M.values
    return ny_df

def compare_order(df,index,use_proba_template=True):
    """
    return the template(int) and 
    indiviual sample score(numpy.array)
    """
    subset=df[:,list(index)]
    ### suppose the index df.iloc[:,0]<df.iloc[:,1]
    ### so if half of df.iloc[:,0]<df.iloc[:,1], we assign 1 ,else 0
    subset1=subset[:,0]-subset[:,1]
    ### count the number smaller than 0
    index=subset1<0
    suppose_rank=np.count_nonzero(index)
    opp_rank=np.count_nonzero(~index)
    if suppose_rank>opp_rank:
        number=1
        ### True:1,False:0
        index=index.astype(int)
    elif suppose_rank<opp_rank:
        number=0
        index=np.logical_not(index).astype(int)
    else:
        number=random.choice([0,1])
        if number==1:
            index=index.astype(int)
        else:
            index=np.logical_not(index).astype(int)
    if use_proba_template:
        if number==1:
            y=Counter(index)
            ### use the most frequently
            number1=y[1]/index.shape[0]
        else:
            y=Counter(index)
            number1=y[0]/index.shape[0]
    #index=pd.DataFrame(index) 
    # return a np.array
    return number1,index

    
def read_target_file(phenotype,cal_control=False,shuffle=False):
    if phenotype=='liver':
        file_name=file_name1
    if phenotype=='Ps':
        file_name=file_name2
    target_df=read_extract(file_name,work_dir,cal_control,shuffle)
    return target_df
def shuffle_total_df(df):
    subset_df=df.loc[:,['label1','phenotype']]
    shulle_res=subset_df.sample(frac=1)
    df.loc[:,['label1','phenotype']]=shulle_res.values
    return df
def get_target_file(phenotype,cal_control):
    total_df11=total_df1.loc[total_df1['phenotype']==phenotype]
    if cal_control:
        total_df11=total_df11.loc[total_df11['label1']==0]
    else:
        total_df11=total_df11.loc[total_df11['label1']==1]
    ###remove first ,last-1,last colunm
    column_indices=[0,len(total_df11.columns)-2,len(total_df11.columns)-1]
    drop_names=total_df11.columns[column_indices]
    total_df_M=total_df11.drop(drop_names,axis=1)
    ny_df=total_df_M.values
    return ny_df
def get_base_score(phenotype,cal_control=False,shuffle=False,total_permut=False):
    """
    return tuple:(score,template)
    """
    if total_permut:
        target_df=get_target_file(phenotype,cal_control)
    else:
        target_df=read_target_file(phenotype,cal_control,shuffle)
    comb_iterable=get_combinations()
    final=calculate_score(comb_iterable,target_df,phenotype)
    return final

def calculate_score(comb_iterable,target_df,phenotype):
    """
    return a int and a np.array
    """
    template=[]
    dimension=target_df.shape[0]
    sample_score=np.empty((dimension,1))
    if phenotype=='liver':
        cache_dict=cache_dict1
    else:
        cache_dict=cache_dict2
    while True:
        try:
            index=next(comb_iterable)
            if index not in cache_dict: 
                number,score=compare_order(target_df,index)
                cache_dict[index]=(number,score)
            else:
                number,score=cache_dict[index]
            template.append(number)
            ### convert into two dimension
            score=np.reshape(score,(score.shape[0],-1))
            sample_score=np.concatenate([sample_score,score],axis=1)
        except StopIteration:
            break
    ### remove the first column
    sample_score=sample_score[:,1:]
    n_combinations=sample_score.shape[1]
    print(n_combinations)
    ### get normalized sample_scores, columns:combinations, row:samples
    sample_score_Normal=sample_score.sum(axis=1)/n_combinations
    final_score=sample_score_Normal.mean(axis=0)
    template=np.array(template)
    return final_score,template

def generate_possible_combinations(forward=False):
    """
    return result:
    {3:{{(2,7,9):iterable}},
    4:{},...22:{}}
    3: the round
    (2,7,9): index
    iterable: combinations that select2
    """
    all_combninations=defaultdict(list)
    ### gene set size from 3 to 22
    if forward:
        range_list=[i for i in range(3,24)]
    else:
        range_list=[i for i in range(23,2,-1)]
    for i in range_list:
        combinations1=get_combinations(r=i)
        small_combinations=defaultdict()
        while True:
            try:
                one_comb=next(combinations1)
                ### one_comb: tuple
                res_inter=get_combinations(one_comb,2)
                small_combinations[one_comb]=res_inter
            except StopIteration:
                break
        all_combninations[i]=small_combinations
    return all_combninations
    

def result_1(cal_control=False,shuffle=False,use_proba_template=True,total_permut=False):
    scores,template1=get_base_score('liver',cal_control,shuffle,total_permut)
    scores1,template2=get_base_score('Ps',cal_control,shuffle,total_permut)
    print(scores,scores1)
    mean_score=(scores+scores1)/2
    if use_proba_template:
        #temp_dist=np.linalg.norm(template1-template2)
        #temp_similar=temp_dist
        cos_sim=template1.dot(template2)/(np.linalg.norm(template1)*np.linalg.norm(template2))
        temp_similar=cos_sim
    else:
        temp_distance=(template1==template2)
        Temp_distance=temp_distance.sum()
        temp_similar=Temp_distance/len(temp_distance)
    return temp_similar,mean_score  






if __name__=='__main__':
    work_dir='D:/two_diseases/'
    test=True
    if test:
        file_name1='LIHC_for_python.csv'
        file_name2='Ps_for_python.csv'
    else:
        file_name1='validation_liver.csv'
        file_name2='validate_ps.csv'
    key_words='total_permutation'
    cache_dict1=dict()
    cache_dict2=dict()

    if key_words=='total_permutation':
        if test:
            total_df=pd.read_csv('D:/two_diseases/permutation_df.csv')
        else:
            total_df=pd.read_csv('D:/two_diseases/permutation_df_validate.csv')
        res=[]
        for _ in range(5000):
            total_df1=shuffle_total_df(total_df)
            cache_dict1=dict()
            cache_dict2=dict()
            similar1,_=result_1(cal_control=False,shuffle=True,total_permut=True)
            ### relase the cache_dict
            cache_dict1=dict()
            cache_dict2=dict()
            similar2,_=result_1(cal_control=True,shuffle=True,total_permut=True)
            result=abs(similar1-similar2)
            print(similar1,similar2)
            res.append(result)
        df=pd.DataFrame(res)
        
        if test:
            save_dir=os.path.join(work_dir,'permutation_test_all.csv')
        else:
            save_dir=os.path.join(work_dir,'permutation_test_all_validate.csv')
        df.to_csv(save_dir)
    else:
        similar1,_=result_1()
        cache_dict1=dict()
        cache_dict2=dict()
        similar2,_=result_1(cal_control=True)
        print(similar1,similar2)
        print(abs(similar1-similar2))

    


    
     


