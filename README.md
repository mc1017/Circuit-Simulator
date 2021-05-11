# electronicsdesignproject
Summer Term Electronics Design Project 
Team Members : Marco, Weihan , James

Start Date: 11/5/2021
End Date: 13/6/2021


Github Naming Convention:

Stable - stable	- Accepts merges from Working and Hotfixes
Working - master - Accepts merges from Features/Issues and Hotfixes
Features/Bug topic-* - Always branch off HEAD of Working
Hotfix	hotfix-*	- Always branch off Stable

If you are working on a new feature i.e. Solving inverse of conductance matrix, start a new branch named feature/solveconductancematrix 

Side note:
1. Head is a pointer that points to the latest commit of a branch. Always create a new branch of the HEAD of master when trying to add a new feature

2. DO NOT use hotfix for non-urgent bugs. Hotfix is only used to fix bugs off the Stable version of program. Use branch (bug)

If unsure of the convention, always refer to the link below!
https://gist.github.com/digitaljhelms/4287848
(There are a lot of naming conventions out there. This is the one we are using for this project)


Commit Messages convention:

<type>[optional scope]:<description>
etc. fix: solve inverse matrix
which means the feature of solve inverse matrix has a bug and is fixed

info - used for info and comments
fix - used for bug fixes
feat - used for feature
hotfix - used for bug fixes off stable


