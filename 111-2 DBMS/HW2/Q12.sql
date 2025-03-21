SELECT E.FNAME, E.LNAME 
FROM EMPLOYEE AS E 
WHERE EXISTS 
(SELECT *
FROM DEPENDENT AS D
WHERE E.SSN = D.ESSN AND E.SEX = D.SEX AND E.FNAME = D.DEPENDENT_NAME);