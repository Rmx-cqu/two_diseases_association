import itertools
import pandas as pd
x=[0.1308017,0.1675258, 0.1417526, 0.1082474,
 0.1365979,0.1108247, 0.2800000, 0.5000000,
0.4473684,0.3421053,0.4210526,0.3684211]
Res=[]
x1=sum(x[0:6])/6
x2=sum(x[6:12])/6
base_line=abs(x1-x2)
print(base_line)

for i_tuple in itertools.combinations(x,6):
    y=[ i for i in x  if i not in i_tuple ]
    first=list(i_tuple)
    first_mean=sum(first)/len(first)
    second_mean=sum(y)/len(y)
    statistic=abs(first_mean-second_mean)
    Res.append(statistic)

count=0
total=len(Res)
save_dict={'test':Res}
df=pd.DataFrame(save_dict)
df.to_csv('D:/two_diseases/pathways_permutate.csv')
for i in Res:
    if i>base_line:
        count+=1
print(count)
print((count+1)/(total+1))  
