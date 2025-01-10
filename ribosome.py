import re
from functools import reduce


def read_codons(codon_file):
  # This will open a file with a given path (codon_file)
  file = open(codon_file)
  global myDic
  myDic = {}
  global myRange
  myRange = []
  global mx
  mx = 0
  global mn
  mn = 0


  # Iterates through a file, storing each line in the line variable
  for line in file:
    myValues = line.split(": ")[0]
    myKey = line.split(": ")[1]

    #if there is only one sequance for the value
    myKey = patternMatch(myKey.replace(" ", "").replace("\n", ""))
    if line.__contains__(',') != True:
      if myKey in myDic:
        continue
      myDic[myKey] = myValues
    else:
      myKey = myKey.split(',')
      for i in myKey:
        if i in myDic:
          continue
        myDic[i.replace(" ", "").replace("\n", "")] = myValues

  #remove any incorrect value
  checkNum()
  #remove repitetion
  myRange = list(set(myRange))
  mx = max(myRange)
  mn = min(myRange)


def read_evals(eval_file):
  # This will open a file with a given path (eval_file)
  file = open(eval_file)
  global myEval
  myEval = []
  reading = []
  represent = []

  # Iterates through a file, storing each line in the line variable
  for line in file:
    # Insert code here
    values = line.split(": ")[1]
    reading = values.split(", ")[0]
    represent = values.split(", ")[1]
    myEval.append([reading, represent.replace("\n", "")])


def encode(sequence):
  if sequence in myDic.values():
    for myKey in myDic:
      if myKey == sequence:
        return myKey
  aminoAcid = sequence.split(" ")
  c = ""
  value = ""
  for amino in aminoAcid:
    for myKey in myDic:
      if myDic[myKey] == amino:
        value = myKey
        c += str(value)
        break
  return str(c)


def decode(sequence):
  aminoAcids = ""  #to return
  i = 0
  
  if sequence in myDic:
    return myDic[sequence]

  if len(myRange) == 1:
    while i+3 <= len(sequence):
      myRNA = sequence[i:i+3]
      if myRNA in myDic:
        aminoAcids += myDic[myRNA] + " "
        i+=3
      else:
        i+=1

  else:
    while i < len(sequence):
      x = False
      for j in reversed(range(int(mn), int(mx + 1))):
        myRNA = sequence[i:i + j]
        if myRNA in myDic:
          aminoAcids += myDic[myRNA] + " "
          i += j - 1
          x = True
          break
      if not x:
        i += 1
  return aminoAcids[:-1]

def operate(sequence, eval_name):
  #variables----------------------------------------
  op_num = int(eval_name[len(eval_name) - 1]) - 1
  op = myEval[op_num][1]
  opE = myEval[op_num][0]
  newDecode_RNA = []
  startBuilding = False
  returnRNA = ""
  #-------------------------------------------------

  #which way to go?
  #print(f"{eval_name} = {op} + {opE} ")
  if opE == "R":
    #reverse string
    RNA = sequence[len(sequence)::-1]

  else:
    RNA = sequence

  #the decoded list to code from
  decode_RNA = list(decode(RNA).split(" "))

  #start operating
  #Prefix notation (+ 1 2)
  #---------------------------------------------------------------
  if op == "PR":
    i = -1
    while i < len(decode_RNA) - 1:
      i += 1
      item = decode_RNA[i]
      #print(f"[{startBuilding}]>{item}<")

      if item == "START":
        startBuilding = True
        continue

      if item == "STOP":
        startBuilding = False
        continue

      if startBuilding and item == "EXCHANGE":
        newDecode_RNA.append("=" + decode_RNA[i + 1])
        i += 1
        continue

      if startBuilding and item == "SWAP":
        newDecode_RNA.append(decode_RNA[i + 2])
        newDecode_RNA.append(decode_RNA[i + 1])
        i += 2
        continue

      if startBuilding and item == "DEL":
        i += 1
        continue

      if startBuilding:
        newDecode_RNA.append(item)
        #print(f"=>{newDecode_RNA[len(newDecode_RNA)-1]}<=")
    #end of for loop
    #encode
    for amino in newDecode_RNA:
      y = False
      if "=" in amino:
        for key, value in myDic.items():
          if value == amino[1:] and y:
            returnRNA += key
            continue
          if value == amino[1:]:
            y = True
            continue
      else:
        returnRNA += encode(amino)
    return returnRNA
    #print(f"{RNA} \n\n= {decode_RNA} \n\n=> {newDecode_RNA} \n\n=> {returnRNA}")

  #---------------------------------------------------------------

  #Postfix notation (1 2 +)
  if op == "PO":
    i = -1
    while i < len(decode_RNA) - 1:
      i += 1
      item = decode_RNA[i]
      #print(f"[{startBuilding}]>{item}<")

      if item == "START":
        startBuilding = True
        continue

      if item == "STOP":
        startBuilding = False
        continue

      if startBuilding and item == "EXCHANGE":
        s = "=" + newDecode_RNA[len(newDecode_RNA) - 1]
        newDecode_RNA[len(newDecode_RNA) - 1] = s
        #i += 1
        continue

      if startBuilding and item == "SWAP":
        if len(newDecode_RNA) >= 2:
          p1 = newDecode_RNA.pop()
          p2 = newDecode_RNA.pop()
          newDecode_RNA.append(p1)
          newDecode_RNA.append(p2)
        continue

      if startBuilding and item == "DEL":
        if len(newDecode_RNA) != 0:
          newDecode_RNA.pop()
        continue

      if startBuilding:
        newDecode_RNA.append(item)

    #end of for loop

    #encode
    for amino in newDecode_RNA:
      y = False
      if "=" in amino:
        for key, value in myDic.items():
          if value == amino[1:] and y:
            returnRNA += key
            continue
          if value == amino[1:]:
            y = True
            continue
      else:
        returnRNA += encode(amino)

    #print(f"{RNA} \n\n= {decode_RNA} \n\n=> {newDecode_RNA} \n\n=> {returnRNA}")
    return returnRNA

  #---------------------------------------------------------------

  #Infix notation (1 + 2)
  if op == "I":
    i = -1
    while i < len(decode_RNA) - 1:
      i += 1
      item = decode_RNA[i]
      #print(f"[{startBuilding}]>{item}<")

      if item == "START":
        startBuilding = True
        continue

      if item == "STOP":
        startBuilding = False
        continue

      if startBuilding and item == "EXCHANGE":
        # s = "=" + newDecode_RNA[len(newDecode_RNA)-1]
        # newDecode_RNA[len(newDecode_RNA)-1] = s
        # #i += 1
        # continue
        newDecode_RNA.append("=" + decode_RNA[i + 1])
        i += 1
        continue

      if startBuilding and item == "SWAP":
        temp = newDecode_RNA.pop()
        newDecode_RNA.append(decode_RNA[i + 1])
        newDecode_RNA.append(temp)
        i += 1
        continue

      if startBuilding and item == "DEL":
        i += 1
        continue

      if startBuilding:
        #print(f"=>{item}<=")
        newDecode_RNA.append(item)

    #end of for loop

    #encode
    for amino in newDecode_RNA:
      y = False
      if "=" in amino:
        for key, value in myDic.items():
          if value == amino[1:] and y:
            returnRNA += key
            continue
          if value == amino[1:]:
            y = True
            continue
      else:
        returnRNA += encode(amino)

    #print(f"{RNA} \n\n= {decode_RNA} \n\n=> {newDecode_RNA} \n\n=> {returnRNA}")
    return returnRNA

#My methods------------------------------------------------------


def pattern(match):
  content = match.group(1)
  num = int(content)
  return match.group(0)[0] * num


def patternMatch(str):
  return re.sub(r'[A-Z]{(\d+)}', pattern, str)


#removes any code with a number on it
def checkNum():
  toCheck = re.compile(r"[0-9]")
  newDict = {}

  for myKey in myDic:
    matched_key = toCheck.search(myKey)
    matched_value = toCheck.search(myDic[myKey])
    filter_key = len(list(filter(lambda c: c in "ACGU", myKey))) == len(myKey)
    if not matched_key and not matched_value and filter_key:
      myRange.append(len(myKey))
      newDict[myKey] = myDic[myKey]

  myDic.clear()
  myDic.update(newDict)