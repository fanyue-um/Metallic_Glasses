from random import shuffle

numStruct = 1000
train_number = 0
validation_number = 0
testing_number = 1000
parameter_number = 0

for inilabel in ['11_1', '12_1']:
    if inilabel == '12_1' or inilabel == '11_1':
        sublabelSet = ['']

    for sublabel in sublabelSet:
        coolingrate = inilabel + sublabel

        structures = list(range(numStruct))
        dataname = './SOAP_train_' + coolingrate + '.txt'
        with open(dataname,'r') as data:
            lines = data.readlines()        

        shuffle(structures)
        with open('../SOAP_training.txt','a') as training:
            with open('../ML_training_ISenergy.txt','a') as trainingenergy:
                for i in range(train_number):
                    training.write('%s' % lines[structures[i]])
                    line = lines[structures[i]].split()
                    trainingenergy.write('%s\n' % line[0])

        with open('../SOAP_validation.txt','a') as validation:
            with open('../ML_validation_ISenergy.txt','a') as validationenergy:
                for i in range(train_number, train_number + validation_number):
                    validation.write('%s' % lines[structures[i]])
                    line = lines[structures[i]].split()
                    validationenergy.write('%s\n' % line[0])

        with open('../SOAP_testing.txt','a') as testing:
            with open('../ML_testing_ISenergy.txt','a') as testingenergy:
                for i in range(train_number + validation_number, train_number + validation_number + testing_number):
                    testing.write('%s' % lines[structures[i]])
                    line = lines[structures[i]].split()
                    testingenergy.write('%s\n' % line[0])

        structures = structures[:train_number]
        shuffle(structures)
        with open('../SOAP_parameter.txt','a') as parameter:
            with open('../ML_parameter_ISenergy.txt','a') as parameterenergy:
                for i in range(parameter_number):
                    parameter.write('%s' % lines[structures[i]])
                    line = lines[structures[i]].split()
                    parameterenergy.write('%s\n' % line[0])
