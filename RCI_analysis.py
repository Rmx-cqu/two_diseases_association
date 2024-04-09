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
    
def get_maximize_scores():
    """
    return a dataframe recorded the process and a optimium result
    """
    df=read_target_file(phenotype='liver',cal_control=False)
    df1=read_target_file(phenotype='Ps',cal_control=False)
    all_comb=generate_possible_combinations()
    ### set the default values
    max_temp_simiar=0
    minum_distance=1
    previous_mean=0
    previous_opti=0
    colnames=['N_genes','gene_index','distance','similar','rank_conserv_1','rank_conserv_2','opti_labmda']
    info_df=pd.DataFrame(columns=colnames)
    for key,value in all_comb.items():
        ### key : N_genes(round)
        ### recorded the local optimial of each round
        local_minum_dist=1
        local_max_simiar=0
        local_mean=0
        local_previous_opti=0
        # local_key=0
        # save_final_score=0
        # save_final_score1=0
        # save_opti_lambda=0
        for key1,value1 in value.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')
            final_score1,template1=calculate_score(value2,df1,'Ps')
            
            ### calculate the distance
            distance=abs(final_score-final_score1)

            # ### calculate the template similar
            # temp_distance=(template==template1)
            
            # ### count the number of true
            # Temp_distance=temp_distance.sum()
            
            # ### normalize the count
            # temp_similar=Temp_distance/len(temp_distance)
            
            cos_sim=template1.dot(template)/(np.linalg.norm(template1)*np.linalg.norm(template))
            temp_similar=cos_sim
            ### save the result through the row direction
            # info=[key,key1,distance,temp_similar]
            # info=pd.DataFrame(info)
            # info=info.T
            # info.columns=colnames
            # info_df=pd.concat([info_df,info],axis=0)
            
            ### how to fresh the result
            mean=(final_score+final_score)/2
            paramter=0.6
            ### defince a optimal vauel
            optimize_value=paramter*mean+(1-paramter)*temp_similar
            threshold=0.005
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                
                max_temp_simiar=temp_similar
                minum_distance=distance
                optimum_result=(select_index,minum_distance,
                    max_temp_simiar)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value
                
                local_max_simiar=temp_similar
                local_minum_dist=distance
                local_key=key1
                #local_mean=mean
                save_final_score=final_score
                save_final_score1=final_score1
                save_opti_lambda=local_previous_opti

            
        ### save the result through the row direction
        info=[key,local_key,local_minum_dist,local_max_simiar,save_final_score,save_final_score1,save_opti_lambda]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0)         
    return info_df,optimum_result
def search_sub_structure_one_phenotype():
    """
    return a dataframe recorded the process and a optimium result
    """
    """
    test in TCGA-LIHC
    """
    df=read_target_file(phenotype='liver',cal_control=False)
    #df1=read_target_file(phenotype='Ps',cal_control=False)
    all_comb=generate_possible_combinations()
    ### set the default values
    previous_opti=0
    colnames=['N_genes','gene_index','rank_conserv','opti_value']
    info_df=pd.DataFrame(columns=colnames)
    for key,value in all_comb.items():
        ### key : N_genes(round)
        ### recorded the local optimial of each round

        local_previous_opti=0
        for key1,value1 in value.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')
            
            ### how to fresh the result
            ### defince a optimal vauel
            optimize_value=final_score
            threshold=0.001
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                optimum_result=(select_index,previous_opti)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value
                
                local_key=key1
                #local_mean=mean
                save_final_score=final_score
                save_opti_val=local_previous_opti     
        ### save the result through the row direction
        info=[key,local_key,save_final_score,save_opti_val]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0)         
    return info_df,optimum_result

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

def comnibations_with_fixed_group(fixed_group,remaing_elements):
    for i in remaing_elements:
        ###generate a tuple
        combine_group=tuple(fixed_group)+(i,)
        yield combine_group

def forward_selection():
    ### forward_select
    df=read_target_file(phenotype='liver',cal_control=False)
    df1=read_target_file(phenotype='Ps',cal_control=False)
    start_length=3
    list_length=23
    all_list=[i for i in range(list_length)]
    max_temp_simiar=0
    minum_distance=1
    previous_mean=0
    previous_opti=0
    colnames=['N_genes','gene_index','distance','similar','rank_conserv_1','rank_conserv_2','opti_labmda']
    info_df=pd.DataFrame(columns=colnames)
    for round in range(start_length,list_length):
        numb_genes=round
        small_comb=defaultdict()
        if round==start_length:
            combinations1=get_combinations(all_list,r=numb_genes)
        else:
            remaing_elements=[x for x in all_list if x not in choose_list ]
            combinations1=comnibations_with_fixed_group(choose_list,remaing_elements)
        while True:
            try:
                one_comb=next(combinations1)
                ### one_comb: tuple
                res_inter=get_combinations(one_comb,2)
                small_comb[one_comb]=res_inter
            except StopIteration:
                break
        ### just run one time
        local_minum_dist=1
        local_max_simiar=0
        local_mean=0
        local_previous_opti=0
        for key1,value1 in small_comb.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')
            final_score1,template1=calculate_score(value2,df1,'Ps')
            
            ### calculate the distance
            distance=abs(final_score-final_score1)

            ### calculate the template similar
            temp_distance=(template==template1)
            
            ### count the number of true
            Temp_distance=temp_distance.sum()
            
            ### normalize the count
            temp_similar=Temp_distance/len(temp_distance)

            mean=(final_score+final_score)/2
            paramter=0.5
            ### defince a optimal vauel
            optimize_value=paramter*mean+(1-paramter)*temp_similar
            threshold=0.05
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                max_temp_simiar=temp_similar
                minum_distance=distance
                optimum_result=(select_index,minum_distance,
                    max_temp_simiar)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value
                local_max_simiar=temp_similar
                local_minum_dist=distance
                local_key=key1
                #local_mean=mean
                save_final_score=final_score
                save_final_score1=final_score1
                save_opti_lambda=local_previous_opti
        ### save the result through the row direction
        info=[numb_genes,local_key,local_minum_dist,local_max_simiar,save_final_score,save_final_score1,save_opti_lambda]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0) 
        choose_list=list(local_key) 
    return info_df,optimum_result
def forward_selection_one_cohort():
    ### forward_select
    df=read_target_file(phenotype='liver',cal_control=False)
    
    start_length=3
    list_length=23
    all_list=[i for i in range(list_length)]

    previous_mean=0
    previous_opti=0
    colnames=['N_genes','gene_index','rank_conserv','opti_value']
    info_df=pd.DataFrame(columns=colnames)
    for round in range(start_length,list_length):
        numb_genes=round
        small_comb=defaultdict()
        if round==start_length:
            combinations1=get_combinations(all_list,r=numb_genes)
        else:
            remaing_elements=[x for x in all_list if x not in choose_list ]
            combinations1=comnibations_with_fixed_group(choose_list,remaing_elements)
        while True:
            try:
                one_comb=next(combinations1)
                ### one_comb: tuple
                res_inter=get_combinations(one_comb,2)
                small_comb[one_comb]=res_inter
            except StopIteration:
                break
        ### just run one time

        local_mean=0
        local_previous_opti=0
        for key1,value1 in small_comb.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')
            

            ### defince a optimal vauel
            optimize_value=final_score
            threshold=0.05
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                optimum_result=(select_index,previous_opti)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value

                local_key=key1
                #local_mean=mean
                save_final_score=final_score

                save_opti_value=local_previous_opti
        ### save the result through the row direction
        info=[numb_genes,local_key,save_final_score,save_opti_value]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0) 
        choose_list=list(local_key) 
    return info_df,optimum_result


def backward_elimation(list_length=23):
    ### just get the current optimal(will minimize the costs)
    df=read_target_file(phenotype='liver',cal_control=False)
    df1=read_target_file(phenotype='Ps',cal_control=False)
    ### set the default values
    drop_list=[i for i in range(list_length)]
    max_temp_simiar=0
    minum_distance=1
    previous_mean=0
    previous_opti=0
    colnames=['N_genes','gene_index','distance','similar','rank_conserv_1','rank_conserv_2','opti_labmda']
    info_df=pd.DataFrame(columns=colnames)
    for round in range(list_length,3,-1):
        numb_genes=round-1
        combinations1=get_combinations(drop_list,r=numb_genes)
        small_comb=defaultdict()
        while True:
            try:
                one_comb=next(combinations1)
                ### one_comb: tuple
                res_inter=get_combinations(one_comb,2)
                small_comb[one_comb]=res_inter
            except StopIteration:
                break
        ### just run one time
        local_minum_dist=1
        local_max_simiar=0
        local_mean=0
        local_previous_opti=0
        for key1,value1 in small_comb.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')
            final_score1,template1=calculate_score(value2,df1,'Ps')
            
            ### calculate the distance
            distance=abs(final_score-final_score1)

            ### calculate the template similar
            # temp_distance=(template==template1)
            
            # ### count the number of true
            # Temp_distance=temp_distance.sum()
            
            # ### normalize the count
            # temp_similar=Temp_distance/len(temp_distance)
            cos_sim=template1.dot(template)/(np.linalg.norm(template1)*np.linalg.norm(template))
            temp_similar=cos_sim
            mean=(final_score+final_score)/2
            paramter=0.5
            ### defince a optimal vauel
            optimize_value=paramter*mean+(1-paramter)*temp_similar
            threshold=0.1
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                max_temp_simiar=temp_similar
                minum_distance=distance
                optimum_result=(select_index,minum_distance,
                    max_temp_simiar)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value
                local_max_simiar=temp_similar
                local_minum_dist=distance
                local_key=key1
                #local_mean=mean
                save_final_score=final_score
                save_final_score1=final_score1
                save_opti_lambda=local_previous_opti
        ### save the result through the row direction
        info=[numb_genes,local_key,local_minum_dist,local_max_simiar,save_final_score,save_final_score1,save_opti_lambda]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0) 
        drop_list=list(local_key) 
    return info_df,optimum_result

def backward_selection_one_cohort(list_length=23):
    ### just get the current optimal(will minimize the costs)
    df=read_target_file(phenotype='liver',cal_control=False)
    ### set the default values
    drop_list=[i for i in range(list_length)]

    previous_mean=0
    previous_opti=0
    colnames=['N_genes','gene_index','rank_conserv','opti_value']
    info_df=pd.DataFrame(columns=colnames)
    for round in range(list_length,3,-1):
        numb_genes=round-1
        combinations1=get_combinations(drop_list,r=numb_genes)
        small_comb=defaultdict()
        while True:
            try:
                one_comb=next(combinations1)
                ### one_comb: tuple
                res_inter=get_combinations(one_comb,2)
                small_comb[one_comb]=res_inter
            except StopIteration:
                break
        ### just run one time

        local_mean=0
        local_previous_opti=0
        for key1,value1 in small_comb.items():
            select_index=key1
            value2=deepcopy(value1)
            ### values1 : iterable, key1: select index
            final_score,template=calculate_score(value1,df,'liver')

            
            # ### normalize the count

            ### defince a optimal vauel
            optimize_value=final_score
            threshold=0.1
            if (optimize_value-previous_opti)>=threshold:
                previous_opti=optimize_value
                optimum_result=(select_index,previous_opti)
                #previous_mean=mean
            if optimize_value>=local_previous_opti:
                local_previous_opti=optimize_value

                local_key=key1
                #local_mean=mean
                save_final_score=final_score

                save_opti_lambda=local_previous_opti
        ### save the result through the row direction
        info=[numb_genes,local_key,save_final_score,save_opti_lambda]
        info=pd.DataFrame(info)
        info=info.T
        info.columns=colnames
        info_df=pd.concat([info_df,info],axis=0) 
        drop_list=list(local_key) 
    return info_df,optimum_result


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
    if key_words=='search_sub_structure':     
        ## the globa varaible in the form of list,dict can be changed in the function
        cache_dict1=dict()
        cache_dict2=dict()
        df,opti=get_maximize_scores()
        print(opti)
        save_dir=os.path.join(work_dir,'rank_conver_test1_exhaust.csv')
        df.to_csv(save_dir)
    elif key_words=='search_sub_in_one_phenotype':
        cache_dict1=dict()
        search_type='backward'
        if search_type=='exhuast':
            df,opti=search_sub_structure_one_phenotype()
            print(opti)
            save_dir=os.path.join(work_dir,'rank_conver_LIHC_1_27.csv')
        elif search_type=='forward':
            df,opti=forward_selection_one_cohort()
            print(opti)
            save_dir=os.path.join(work_dir,'rank_conver_LIHC_1_27_forward.csv')
        else:
            df,opti=backward_selection_one_cohort()
            print(opti)
            save_dir=os.path.join(work_dir,'rank_conver_LIHC_1_27_backward.csv')
        df.to_csv(save_dir)
    elif key_words=='permutation':
        #### permutation test
        res=[]
        for _ in range(5000):
            cache_dict1=dict()
            cache_dict2=dict()
            similar1,_=result_1(cal_control=False,shuffle=True)
            
            ### relase the cache_dict
            cache_dict1=dict()
            cache_dict2=dict()
            similar2,_=result_1(cal_control=True,shuffle=True)
            result=abs(similar1-similar2)
            print(similar1,similar2)
            res.append(result)
        df=pd.DataFrame(res)
        save_dir=os.path.join(work_dir,'permutation_validate.csv')
        df.to_csv(save_dir)
        print(np.quantile(res,0.999))
    elif key_words=='total_permutation':
        validate=True
        if validate:
            total_df=pd.read_csv('D:/two_diseases/permutation_df_validate.csv')
        else:
            total_df=pd.read_csv('D:/two_diseases/permutation_df.csv')
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
        
        if validate:
            save_dir=os.path.join(work_dir,'permutation_test_all_validate.csv')
        else:
            save_dir=os.path.join(work_dir,'permutation_test_all.csv')
        df.to_csv(save_dir)

    elif key_words=='back_elimination':
        df,opti=backward_elimation()
        print(opti)
        save_dir=os.path.join(work_dir,'rank_conver_back.csv')
        df.to_csv(save_dir)
    elif key_words=='forward_selection':
        df,opti=forward_selection()
        print(opti)
        save_dir=os.path.join(work_dir,'rank_conver_forward.csv')
        df.to_csv(save_dir)
    else:
        similar1,_=result_1()
        cache_dict1=dict()
        cache_dict2=dict()
        similar2,_=result_1(cal_control=True)
        print(similar1,similar2)
        print(abs(similar1-similar2))

    


    
     


