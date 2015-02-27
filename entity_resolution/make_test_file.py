f0 = open('er-test.db')
f1 = open('er-test-eclipse.db','w')
class_dic = {}
class_count = 0
author_dic = {}
author_count = 0
venue_dic = {}
venue_count = 0
title_dic = {}
title_count = 0
word_dic = {}
word_count = 0
linecount = 0
for line in f0:
	linecount += 1
	print linecount
	if line.startswith('SameBib'):
		break
	line = line.rstrip()
	if not line:
		continue
	[predName,rest] = line.split('(')
	line = rest
	[arguments,rest] = line.split(')')
	fields = arguments.split(',')
	f1.write(predName+"(")
	for field_index in range(len(fields)):
		field = fields[field_index]
		field = field.strip()
		if field.startswith('Class'):
			if not field in class_dic:
				class_dic[field] = class_count
				class_count += 1
			f1.write(str(class_dic[field]))	
		elif field.startswith('Author'):
			if not field in author_dic:
				author_dic[field] = author_count
				author_count += 1
			f1.write(str(author_dic[field]))	
		elif field.startswith('Venue'):
			if not field in venue_dic:
				venue_dic[field] = venue_count
				venue_count += 1
			f1.write(str(venue_dic[field]))	
		elif field.startswith('Title'):
			if not field in title_dic:
				title_dic[field] = title_count
				title_count += 1
			f1.write(str(title_dic[field]))	
		elif field.startswith('Word'):
			if not field in word_dic:
				word_dic[field] = word_count
				word_count += 1
			f1.write(str(word_dic[field]))	
		if field_index != len(fields)-1:
			f1.write(",")
	f1.write(")\n")
	#print class_dic
	#print author_dic
	#print venue_dic
	#print title_dic
	#print word_dic
	#input()
print class_count
print author_count
print venue_count
print title_count
print word_count
f0.close()
f1.close()
