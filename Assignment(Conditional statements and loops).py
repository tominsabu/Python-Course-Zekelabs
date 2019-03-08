# 1.	Write a Python program to find those numbers which
#       are divisible by 7 and multiple of 5, between 1500 and 2700 (both included). 

# count=0
# for value in range(1500,2701):
# 	if value%7==0 and value%5==0:
# 		print(value,end=' ')
# 		count+=1
# print('number of values=',count)	


# 2.	Write a Python program to convert temperatures to and from celcius, fahrenheit 
#       [ Formula : c/5 = f-32/9 [ where c = temperature in celsius and f = temperature in fahrenheit ] 
    

# f=100
# c=5*(f-(32))/9
# print('The value of temperature in Celsius is',c)
# # s='The value of temperature in Celsius is {}'.format(c)
# # print(s)

# c=37.77777777777778
# f=32+(9/5)*c
# print('The value of temperature in Fahrenheit is',f)
# # s1='The value of temperature in Fahrenheit is {}'.format(f)
# # print(s1)


# 3.	Write a Python program to count the number of even and odd numbers from a series of numbers.  
# Sample numbers : numbers = (1, 2, 3, 4, 5, 6, 7, 8, 9) 
# Expected Output : 
# Number of even numbers : 5
# Number of odd numbers : 4 

# numbers = (1, 2, 3, 4, 5, 6, 7, 8, 9)
# count_even=0
# count_odd=0
# for value in numbers:
# 	if value%2==0:
# 		count_even+=1
# 	else:
# 		count_odd+=1
# print('Number of even numbers =',count_even,',and the number of odd numbers =',count_odd)



# 4.	How many times will the following code print "Welcome to Python"?

# count = 0
# while count < 10:
#     print("Welcome to Python")
#     count += 1

# Ans: It counts from 0 to 9, a total number of 10 times


# 5.	What is the output of the following code?

# x = 0
# while x < 4:
#     x = x + 1

# print("x is", x)

# Ans:4
# x    : 0 1 2 3 4
# x=x+1: 1 2 3 4


# 6.	What will be displayed when the following code is executed?

# number = 6
# while number > 0:
#     number -= 3
#     print(number)

# Ans: 3,0

# number         :  6 3
# number=number-3:  3 0


# 7.	Write a Python program which iterates the integers from 1 to 50. 
# For multiples of three print "Fizz" instead of the number and for the multiples of five print "Buzz". 
# For numbers which are multiples of both three and five print "FizzBuzz".
# Sample Output : 
# fizzbuzz
# 1
# 2
# fizz
# 4 
# buzz 
# for value in range(1,51):
# 	if value%3==0 and value%5==0:
# 		print('FizzBuzz')
# 	elif value%3==0:
# 		print('Fizz')
# 	elif value%5==0:
# 		print('Buzz')
# 	else:
# 		print(value)


# 8.	Write a Python program to calculate the sum and average of n integer numbers.
	
# n=100
# sum=0
# count=1
# while count<=n:
# 	sum=sum+count
# 	count+=1
# print('sum=',sum)
# print('Average=',sum/n)


# 9.	Write a Python program to get the Fibonacci series between 0 to 50. 

# n1=0
# n2=1
# while n1<=50:
# 	print(n1)
	# n2,n1=n1+n2,n2


# 10.	What will be the output for :
x = 'abcd'
for i in range(len(x)):
    x = 'a'
    print(x)

# Ans: four as	
# i: 1 2 3 4
# x: a a a a

	
