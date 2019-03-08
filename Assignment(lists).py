# 1. Write a Python function called compare date that takes as arguments two lists of two integers each. 
# Each list contains a month and a year, in that order. 
# The function should return -1 if the first month and year are earlier than the second month and year, 
# 0 if they are the same, and 1 if the first month and year are later than the second. 
# Your code should work for any legal input for month and year. 
# Example calls and expected output are shown below:
# >>> compare_date( [10,1995], [8,1995] )
# 1
# >>> compare_date( [5,2010], [5,2010] )
# 0
# >>> compare_date( [10,1993], [8,1998] )
# -1

#Answer:

# def compare_date(month1,year1,month2,year2):
# 	l1=[month1,year1]
# 	l2=[month2,year2]
# 	print(l1,l2)
# 	if (l1[0]<=l2[0] and l1[1]<l2[1]) or (l1[0]<l2[0] and l1[1]==l2[1]):
# 		print(-1)
# 	elif l1[0]==l2[0] and l1[1]==l2[1]:
# 		print(0)
# 	else:
# 		print(1)
# compare_date(12,1995,10,1995)

#--------------------------------------------------------------------------------------------------------------



# 2. Assume v is a list containing numbers. 
# Write Python code to find and print the highest two values in v. 
# If the list contains only one number, print only that number. 
# If the list is empty, print nothing. For example, if we assigned 
# v = [ 7, 3, 1, 5, 10, 6 ]
# then the output of your code should be something like
# 7 10
# If we are given that
# v = [ 7 ]
# then the output of your code should be 7

#Answer

# v = [4,5,4,44,6,6,6,12,56]
# max=0
# sec_max=0
# if len(v)==1:
# 	print(v)
# elif len(v)>1 and v==v[::-1]:     
# 	print(v)                       
# elif len(v)>1:                    
# 	for value in v:
# 		if value > max:
# 			sec_max=max
# 			max=value
# 		elif value > sec_max and value < max:
# 			sec_max=value
# 	print(sec_max,max)

#-----------------------------------------------------------------------------------------------------------------


# 3. Consider a data, where just the name of the restaurant, the type of restaurant, and the ratings are provided. 
# Assume these values have already been read into a list of lists of the form below:
# restaurants = [ [ Acme, Italian, 2, 4, 3, 5],[ Flintstone, Steak, 5, 2, 4, 3, 3, 4],[ Bella Troy, Italian, 1, 4, 5] ]
# Write a segment of Python code that prints all Italian restaurants 
# in the restaurants list that have no ratings of value 1 and at least one rating of value 5. 
# In the above example, Acme would be printed in the output, but Flintstone and Bella Troy would not. 
# Flintstone is not Italian and Bella Troy has a 1 rating. 

#Answer

# restaurants=[ [ 'Acme', 'Italian', 2, 4, 3,5],[ 'Flintstone', 'Steak', 5, 2, 4, 3, 3, 4],[ 'Bella Troy', 'Italian',1, 4, 5]] 
# count1=0
# while count1<len(restaurants):
# 	if restaurants[count1][1]=='Italian':
# 		if (not(1 in restaurants[count1][2:])) and (5 in restaurants[count1][2:]):
# 			print(restaurants[count1][0])
# 	count1+=1

#----------------------------------------------------------------------------------------------------------------




# 4. In the game of chess you can often estimate how well you are doing by adding the values of the pieces
#  you have captured. The pieces are Pawns, Bishops, Knights, Rooks and Queens. Their values are
# P - (P)awn, value = 1
# B - (B)ishop, value = 3
# K - (K)night, value = 3
# R - (R)ook, value = 5
# Q - (Q)ueen, value = 9
# Write a Python function called chess_score that takes a single string as an argument and returns the
# combined  values represented by the pieces in the string. 
# You may assume that only P, B, K, R and Q appear in the string. 
# You may not use any if statements and you may not use any loops. 
# As an example,
# print chess_score( BQBP )
# should output the value 16 because there are 2 Bishops (3 points each), 1 Queen (9 points each), 
# and 1 Pawn (1 point each).

#Answer

# def chess_score(s):
# 	l=[['P',1],['B',3],['K',3],['R',5],['Q',9]]
# 	sum=0
# 	for value1 in s:
# 		count=0
# 		while count<len(l):
# 			if value1 in l[count][0]:
# 				sum=sum+l[count][1]
# 				break
# 			count+=1
# 	print(s)	
# 	print(sum)
# chess_score('BQBP')

#-----------------------------------------------------------------------------------------------------------------




# 5. You are given a file that contains, on each line of input, three integers separated by commas. 
# Write a Python program that sums all of the first integers, the second integers, 
# and the third integers, outputting the resulting sums all on one line, separated by commas. 
# As a simple example, if the input is 
# 2, 5, 7 
# 3, 6, 10
# 1, 2, -3
# 2, 4, 1
# Then the output should be
# 8, 17, 15

#Answer

# l1= [2,5,7]
# l2= [3,6,10]
# l3= [1,2,-3]
# l4= [2,4,1]
# sum1=l1[0]+l2[0]+l3[0]+l4[0]
# sum2=l1[1]+l2[1]+l3[1]+l4[1]
# sum3=l1[2]+l2[2]+l3[2]+l4[2]
# s='{},{},{}'.format(sum1,sum2,sum3)
# print(s)

#----------------------------------------------------------------------------------------------------------------


# 6. Write Python code to generate the following ranges
# (a) (100; 99; 98; : : : ; 0)
# (b) (55; 53; 51; : : : ;-1)
# (c) (3; 5; 7; 9; : : : ; 29)
# (d) (-95;-90;-85; : : : ; 85; 90)

#Answer

# l1=[]
# for value in range(0,101):
# 	l1.append(value)
# print(l1[::-1])

# l1=[]
# for value in range(-2,101):
# 	l1.append(value)
# print(l1[l1.index(55):l1.index(-2):-2])


# l1=[]
# for value in range(0,101):
# 	l1.append(value)
# print(l1[l1.index(3):l1.index(30):2])

# l1=[]
# for value in range(-100,100):
# 	l1.append(value)
# print(l1[l1.index(-95):l1.index(91):5])

#-----------------------------------------------------------------------------------------------------------------



# 7. Write a while loop to add all of the numbers in a list v until it reaches a negative number 
# or until it reaches the end of the list. Store the sum in the variable result. 
# Your code should work for any version of v containing only numbers. 
# For example, the value of result should be 25 after the loop for both of the following lists:
# v = [ 10, 12, 3, -5, 5, 6 ]
# v = [ 0, 10, 3, 6, 5, 1 ]

#Answer

# v = [ 10, 12, 3, -5, 5, 6 ]
# count=0
# result=0
# while count<len(v):
# 	if v[count] >=0:
# 		result=result+v[count]
# 		count+=1
# 	else:
# 		break
# print(result)

#-----------------------------------------------------------------------------------------------------------------




# 8. Write Python code that takes a list of numbers, v, and outputs the positive values 
# that are in v in increasing order, one value per line. 
# If there are no positive values, then the output should be the string 'None'. 
# You may assume there is at least one value in the list. 
# As an example, v = [ 17, -5, 15, -3, 12, -5, 0, 12, 22, -1 ]
# Then the output of your code should be
# 12
# 12
# 15
# 17
# 22
# As a second example, if
# v = [ -17, -5, -15, -3, -12, -5, 0, -12, -22, -1 ]
# then then output should be just
# None

#Answer

# v = [ 17, -5, 15, -3, 12, -5, 0, 12, 22, -1 ]
# l1=[]
# for value in v:
# 	if value>0:
# 		l1.append(value)
# l1.sort()
# if len(l1)>=1:
# 	for value in l1:
# 		print(value)
# else:
# 	print(None)

#-----------------------------------------------------------------------------------------------------------------



# 9. What is the output of the following operations:

# mylist = [1,4,8,12,6]
# x = mylist.sort()          
# print (x)

# mylist = [1,4,8,12,6]
# slice1 = mylist[2:4]
# slice1[0] = 20
# print (slice1)
# print (mylist)

# Answer: 

# none(sort() does not return a value)
# [20,12]
# [1,4,8,12,6]

#-----------------------------------------------------------------------------------------------------------------




# 10. Write a Python for loop to print out the values from the list v that are positive (0 is NOT a positive
# number).

#Answer

# v = [ 17, -5, 15, -3, 12, -5, 0, 12, 22, -1 ]
# for value in v:
# 	if value>0:
# 		print(value)

#-----------------------------------------------------------------------------------------------------------------




# 11. What is the output of the following program?

# def spam(a1,b1,a2,b2):
# 	if (a1 == a2) and (b1 > b2):
# 		return 1
# 	else:
# 		return 0
# def egg(a1,b1,a2,b2):
# 	if (a1 > a2) and (b1 == b2):
# 		return 0
# 	else:
# 		return 1
# a1 = 3
# b1 = 4
# a2 = 6
# b2 = 4
# print(spam(a2, b2, a1, b1))
# print( egg(a1, b1, a2, b2))
# c = spam(a1, b2, a2, b1)
# print (c)
# c += egg(a1, b2, a2, b1)
# print (c)

# Answer

# 0
# 1
# 0
# 1

#-----------------------------------------------------------------------------------------------------------------




# 12.Give the output of each of the following
# (a)


# i = 4
# L = [ 0, 12, 3, 5, 2, -1 ]
# while 0 <= i and i < len(L):          
# 	if L[i] < 0:                       
# 		break
# 	else:
# 		i = L[i]
# print (i, L[i])


# Ans:5 -1

# Explanation
# i    :2,3,5
# L[i] :3,5,-1

#-----------------------------------------

# (b)

# tough = 2
# for i in range(2):
# 	s = 1
# 	for j in range(i, tough):
# 		s += tough
# print (s)
# print (tough)
# tough = s
# print (tough)


# Answer: 

# 3 2 3 (Got)

#-----------------------------------------------------------------------------------------------------------------




# 13. Suppose a list of words in alphabetical order has been assigned to the variable called words. 
# For example,
# we might have the assignment words = [ aardvark , abaka, expedite, experience, shoetrees, tastetest, test ]
# Write code to find and output the first and the last string in words that start 
# and end with the same letter and are at least 8 characters long. 
# You may assume that at least one word in words satisfies this condition. 
# You may write a function if you wish. For the above example, the output should be
# 6
# expedite
# tastetest

# words = [ 'aardvark' , 'abaka', 'expedite', 'experience', 'shoetrees', 'tastetest', 'test' ]
# for value in words:
# 	if len(value)==8 and (value[0]==value[len(value)-1]):
# 		print(value,len(value))

# Answer: 

# expedite 8	

#-----------------------------------------------------------------------------------------------------------------




# 14. Please show the output from the following code?

# def get_min(v):
# 	v.sort()
# 	return v[0]
# def get_max(v):
# 	x = max(v)
# 	return x
# v = [ 14, 19, 4, 5, 12, 8 ]
# if len(v) > 10 and get_min(v) > 6:         
# 	print ("Hello")                         
# 	print (v[0])
# 	print (v[4])
# else:
# 	print ("So long")
# 	print (v[0])
# 	print (v[-1])
# if len(v) < 10 or get_max(v):    
# 	print (get_max(v))
# 	print (v[0])
# 	print (get_min(v))
# 	print (v[0])

#  Answer:

# So long,14,8,19,14,4,4

# ----------------------------------------------------------------------------------------------------------------



# 15. Write code that uses a range (and NO loops) to generate the following lists:
# v0 = [ 10, 9, 8, 7, 6, 5, 4, 3 ]
# v1 = [ -10, -3, 4, 11, 18, 25, 32, 39 ]

# Answer

# l=[]
# for nums in range(3,11):
# 	l.append(nums)
# print(l[::-1])

# l=[]
# for nums in range(-10,40):
# 	l.append(nums)
# print(l[::7])

# ----------------------------------------------------------------------------------------------------------------




# 16. Consider the following list of lists of strings:

# wordy = [ [ 'impala', 'malibu', 'camry', 'jetta'],[ 'zebra', 'impala', 'lion', 'impala', 'malibu', 'zebra' ],[ 'tiger', 'lion', 'cowboy', 'jet', '49er' ],[ ],[ 'five', 'seven', 'nine'] ]

# # (a) Show the output:
# print (wordy[2][1])
# print (wordy[1][2][3])      #Ans:n
# print (len(wordy))
# print (len(wordy[1]))

#Answer

# Got

# -------------------------------------------

# (b) Write a loop to print the last word of each list in wordy, stopping when either an empty list is found
# or when there are no more lists. For the above example, the output should be
# jetta
# zebra
# 49er

# Answer

# count=0
# while count<len(wordy):
# 	if len(wordy[count])!=0:
# 		print(wordy[count][len(wordy[count])-1])
# 	else:
# 		break
# 	count+=1
    
# ----------------------------------------------------------------------------------------------------------------




# 17. Show the output of the following code. Make sure we can determine what is output and what is scratch work.

# def remove_something(z):
# 	z.remove( z[z[0]] )

# v = [ 1, 8, [12, 8], 'hello', 'car' ]
# x = 'salad'
# if len(v[2]) >= 2:
# 	if x > v[3]:       
# 		print ('One')  
# 	if v[0] == 1:      
# 		print ('Three')
# 	else:
# 		print ('Two')
# elif len(v) == 5:      
# 	print ('Six')
# else:
# 	v.append('five')
# 	print ('Ten')		

# remove_something(v)
# print (v[1])
# print (v[2])
# v.append(x)
# print (len(v))
# print(v)

# Answer

# Got

#-----------------------------------------------------------------------------------------------------------------




# 18. You are given in variable x a list of lists represented as an NxN grid in which
#  each list corresponds to one row of the grid. For example, a 4x4 grid is given by:
# x = [[1,2,3,4],[4,3,2,1],[2,1,4,2],[2,1,4,5]]
# Write a piece of code to print the grid in the following format with
#  a vertical and horizontal line right in the middle:
# 1 2 | 3 4
# 4 3 | 2 1
# ----|----
# 2 1 | 4 2
# 2 1 | 4 5

#Answer

# x = [[1,2,3,4],[4,3,2,1],[2,1,4,2],[2,1,4,5]]
# s='''
# {}  {} | {}  {}
# {}  {} | {}  {}
# ---- | ----
# {}  {} | {}  {}
# {}  {} | {}  {}'''.format(x[0][0],x[0][1],x[0][2],x[0][3],x[1][0],x[1][1],x[1][2],x[1][3],x[2][0],x[2][1],x[2][2],x[2][3],x[3][0],x[3][1],x[3][2],x[3][3])
# print(s)

#----------------------------------------------------------------------------------------------------------------



# 19. Suppose you are given the scores of two athletes in various competitions, given in two separate lists. 
# Assume there are unknown number of competitions numbered 1,2,3, etc. 
# and the length of the two lists is the same.
# a1 = [11,8,11,9]
# a2 = [11,9,8,12]
# For example according this to list, both athletes got a score of 11 in competition 1. 
# Print the index of all the competitions in which a2 did better. 
# For example, for the above lists, we would print: a2 is better in 2 4
# If there is no value in which a2 is better, then you should print:
# a2 is never better

# Answer

# a1 = [11,8,11,9]
# a2 = [11,9,18,12]
# count=1
# l1=[]
# l2=[]
# l3=[]
# while count<=len(a1):
# 	if a2[count-1]>a1[count-1]:
# 		l1.append(count)
# 	elif a2[count-1]==a1[count-1]:
# 		l2.append(count)
# 	elif (a2[count-1]<a1[count-1]):
# 		l3.append(count)
# 	count+=1
# if len(l1)!=0:
# 	print('a2 is better in',l1)
# if len(l2)!=0:
# 	print('scores are same in',l2)
# if len(l3)==len(a1):
# 	print('a2 is never better')

# ----------------------------------------------------------------------------------------------------------------




# 20. What is the output from the following code:

# L1 = ['cat', 'dog', 'hawk', 'tiger', 'parrot']
# print (L1[1:-1])        #Ans: ['dog', 'hawk', 'tiger']
# print (L1[1:-2])        #Ans: ['dog', 'hawk']
# print (L1[1:-4])        #Ans: []
# print (L1[1:0])         #Ans: []    contradicts with the direction sense
# print (L1[1:10])        #Ans: ['dog', 'hawk', 'tiger', 'parrot']
# print (L1[::-1])        #Ans: ['parrot', 'tiger', 'hawk', 'dog', 'cat']
# print (L1[1:4:2])       #Ans: ['dog', 'tiger']
# print (L1[::-2])        #Ans: ['parrot', 'hawk', 'cat']

# ----------------------------------------------------------------------------------------------------------------




