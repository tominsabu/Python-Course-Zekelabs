# 1.	What is the output of the following Python program section?
# st='Python Programming' 
# x=st.lower()
# y=x.replace('python','Python')
# print (y)

# Ans: Python programming

#x=python programming
#y=Python programming




# # 2.
# txt='Python is powerful and easy.'
# words=txt.split(' ')
# l=len(words)
# print (l)

# Ans: 5

# split(), splits the string txt when ever it encounters a space
# type(words)=list, it has five elements




# # 3.	What would be the output from the following Python program section?
# st=['P','y','t','h','o','n']
# w=''.join(st)
# print (w)

#correct ans: Python

# w=''.join(st) => join the elements of the list st with no seperation between




# 4.	How to check if a String contains only digits?

# s='12345'
# print(s.isdigit())



# 5.	How to count number of vowels and consonants in a String?

#SEE THE ONE BELOW THE FIRST ONE

# s=[]
# for value in range(ord('a'),ord('z')+1):
# 	if chr(value) not in 'aeiou':
# 		s.append(chr(value))
# # print(s)
# s1=('').join(s)
# # print(s1)



#BETTER LOGIC 

# s='areal dogfight 444$$$ 5684AABC'
# count_vowels=0
# count_consonants=0
# count_upvowels=0
# count_upconsonants=0
# for value in s:
# 	if value in "aeiou" or value in 'AEIOU':
# 		count_vowels+=1
# 		if value in "AEIOU":
# 			count_upvowels+=1
# 	elif value.isalpha():                      #how line 81 is not considering condition in line 77
# 		count_consonants+=1                    #if => and, elif => or
# 		if value.isupper():
# 			count_upconsonants+=1
# print('No. of vowels=',count_vowels,'No. of consonants=',count_consonants)
# print('No. of upcase vowels=',count_upvowels,'No.of upcase consonants=',count_upconsonants)





# 6.	How to replace each given character to other e.g. blank with %20? 

# s='samjoe 395 SAMJOE'
# print(id(s))
# s=s.replace(' ','%20')
# print(s)
# print(id(s))




# 7.	How to check if String is Palindrome?

# s='malayalam'
# if s==s[::-1]:
# 	print('String is Palindrome')
# else:
# 	print('String is not Palindrome')




# 8.	You are given a string. Your task is to capitalize each word of the string.

# s='malayala manoraMA *march2019*' 
# s=s.title()
# print(s)


# 9.	You need to write a function to implement an algorithm which will accept a string of characters 
# and should find the highest occurrence of the character and display it. 
# For example if input is "aaaaaaaaaaaaaaaaabbbbcddddeeeeee" it should return "a".

# def highest_occur_char(s):	
# 	l=[]
# 	max=0
# 	for x in s:
# 		l.append(s.count(x))
# 	for y in l:
# 		if y > max:
# 			max=y
# 	for x in s:
# 		if s.count(x) == max:
# 			hioccur=x
# 	print('The character that occurs the most =',hioccur,'(',max,'times)')

# s='aaaaaaaaaaaaaaaaa,,,,,,,,bbbbcddddeeeeee'
# highest_occur_char(s)



# # 10.	How to count occurrence of a given character in String?

# s1="How to Count ocCurrenCe of a given charaCter in String?"
# print(s1.count('C'))
# print(s1.count('c'))
# print(s1.count(s1))