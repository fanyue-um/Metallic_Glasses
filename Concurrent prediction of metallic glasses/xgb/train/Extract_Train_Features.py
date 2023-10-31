train_number = 6000


with open('../SOAP_training.txt','r') as train:
    with open('./TrainingFeatures.txt','a') as output:
        for _ in range(train_number):
            line = train.readline().split()
            for i in range(1,len(line)):
                feature = line[i].split(':')
                output.write("%s " % feature[1])
            output.write("\n")
