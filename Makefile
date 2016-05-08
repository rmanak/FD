.IGNORE:
SHELL = /bin/bash

SOURCES = about.txt features.txt download.txt install.txt start.txt tutorials.txt references.txt manual.txt 

pagefiles=$(SOURCES:.txt=.html)

default: all

all: $(pagefiles) index.html

index.html: about.html
	/bin/cp -f about.html index.html

footbar.pm: footbar
	./Markdown.pl footbar > footbar.pm

sidebar.pm: sidebar
	./Markdown.pl sidebar > sidebar.pm

%.s: %.txt
	cat $*.txt globalauthor globalkeywords globaldescription > $*.txt2
	./find_code.pl < $*.txt2 > $*.s
	/bin/rm -f *.txt2
	/bin/rm -f tmp_code_123 tmp_code_1234

%.w: %.s
	./Markdown.pl $*.s > $*.w

%.html : %.w template footbar.pm sidebar.pm globalkeywords globalauthor globaldescription
	./preppage.pl template sidebar.pm footbar.pm $*.w | inc_file.pl > $*.html
	grep 'table_of_contents' $*.html && python ./pyhat.py -d. $*.html || echo $*.html generated

	
clean:
	/bin/rm -f $(pagefiles) index.html
	
giveaccess:
	chmod -R a+r *
sync:
	git push origin gh-pages
