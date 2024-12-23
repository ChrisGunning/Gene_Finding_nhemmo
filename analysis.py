import matplotlib.pyplot as plt
import numpy as np

FOLDER = "result_files/"

TRAIN_PREDICT_1 = ["Chromosome1_1.csv", "Chromosome2_1.csv", "Chromosome3_1.csv", 
                   "Chromosome4_1.csv", "Chromosome5_1.csv"]
TRAIN_LABEL_1 = ["Chromosome1_1_real.csv", "Chromosome2_1_real.csv", 
                 "Chromosome3_1_real.csv", "Chromosome4_1_real.csv", 
                 "Chromosome5_1_real.csv"]
TRAIN_PREDICT_3 = ["Chromosome1_3.csv", "Chromosome2_3.csv", "Chromosome3_3.csv", 
                   "Chromosome4_3.csv", "Chromosome5_3.csv"]
TRAIN_LABEL_3 = ["Chromosome1_3_real.csv", "Chromosome2_3_real.csv", 
                 "Chromosome3_3_real.csv", "Chromosome4_3_real.csv", 
                 "Chromosome5_3_real.csv"]
TRAIN_PREDICT_6 = ["Chromosome1_6.csv", "Chromosome2_6.csv", "Chromosome3_6.csv", 
                   "Chromosome4_6.csv", "Chromosome5_6.csv"]
TRAIN_LABEL_6 = ["Chromosome1_6_real.csv", "Chromosome2_6_real.csv", 
                 "Chromosome3_6_real.csv", "Chromosome4_6_real.csv", 
                 "Chromosome5_6_real.csv"]

TEST_PREDICT_1 = ["Chromosome22_1.csv", "Chromosome23_1.csv", "Chromosome24_1.csv"]
TEST_LABEL_1 = ["Chromosome22_1_real.csv", "Chromosome23_1_real.csv", "Chromosome24_1_real.csv"]
TEST_PREDICT_3 = ["Chromosome22_3.csv", "Chromosome23_3.csv", "Chromosome24_3.csv"]
TEST_LABEL_3 = ["Chromosome22_3_real.csv", "Chromosome23_3_real.csv", "Chromosome24_3_real.csv"]
TEST_PREDICT_6 = ["Chromosome22_6.csv", "Chromosome23_6.csv", "Chromosome24_6.csv"]
TEST_LABEL_6 = ["Chromosome22_6_real.csv", "Chromosome23_6_real.csv", "Chromosome24_6_real.csv"]


TRANSFER_PREDICT_1 = ["TransferChromosome1_1.csv", "TransferChromosome2_1.csv", 
                      "TransferChromosome3_1.csv", "TransferChromosome4_1.csv", 
                      "TransferChromosome5_1.csv"]
TRANSFER_LABEL_1 = ["TransferChromosome1_1_real.csv", "TransferChromosome2_1_real.csv", 
                 "TransferChromosome3_1_real.csv", "TransferChromosome4_1_real.csv", 
                 "TransferChromosome5_1_real.csv"]
TRANSFER_PREDICT_3 = ["TransferChromosome1_3.csv", "TransferChromosome2_3.csv", 
                      "TransferChromosome3_3.csv", "TransferChromosome4_3.csv", 
                      "TransferChromosome5_3.csv"]
TRANSFER_LABEL_3 = ["TransferChromosome1_3_real.csv", "TransferChromosome2_3_real.csv", 
                    "TransferChromosome3_3_real.csv", "TransferChromosome4_3_real.csv", 
                    "TransferChromosome5_3_real.csv"]
TRANSFER_PREDICT_6 = ["TransferChromosome1_6.csv", "TransferChromosome2_6.csv", 
                      "TransferChromosome3_6.csv", "TransferChromosome4_6.csv", 
                      "TransferChromosome5_6.csv"]
TRANSFER_LABEL_6 = ["TransferChromosome1_6_real.csv", "TransferChromosome2_6_real.csv",  
                    "TransferChromosome3_6_real.csv", "TransferChromosome4_6_real.csv", 
                    "TransferChromosome5_6_real.csv"]

def file_to_seq(file):
    with open(FOLDER + file) as csvfile:
        for line in csvfile:
            sarr = line.strip().split(',')
            return [int(s) for s in sarr]

def compare_seq_acc(real, pred):
    real_seq = np.array(file_to_seq(real))
    pred_seq = np.array(file_to_seq(pred))[:len(real_seq)]

    correct = real_seq == pred_seq
    accuracy = sum(correct)/len(correct)
    return accuracy

def graph_groups(seq_groups, seq_labels, title, xlabel, ylabel, file):
    accuracies = []
    for seqs in seq_groups:
        accuracies+= [sum(seqs)/len(seqs)]

    plt.figure()
    plt.bar(seq_labels, accuracies)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(file)

def main():
    train_1 = [compare_seq_acc(TRAIN_LABEL_1[i],TRAIN_PREDICT_1[i]) for i in range(5)]
    train_3 = [compare_seq_acc(TRAIN_LABEL_3[i],TRAIN_PREDICT_3[i]) for i in range(5)]
    train_6 = [compare_seq_acc(TRAIN_LABEL_6[i],TRAIN_PREDICT_6[i]) for i in range(5)]

    test_1 = [compare_seq_acc(TEST_LABEL_1[i],TEST_PREDICT_1[i]) for i in range(3)]
    test_3 = [compare_seq_acc(TEST_LABEL_3[i],TEST_PREDICT_3[i]) for i in range(3)]
    test_6 = [compare_seq_acc(TEST_LABEL_6[i],TEST_PREDICT_6[i]) for i in range(3)]

    transfer_1 = [compare_seq_acc(TRANSFER_LABEL_1[i],TRANSFER_PREDICT_1[i]) for i in range(5)]
    transfer_3 = [compare_seq_acc(TRANSFER_LABEL_3[i],TRANSFER_PREDICT_3[i]) for i in range(5)]
    transfer_6 = [compare_seq_acc(TRANSFER_LABEL_6[i],TRANSFER_PREDICT_6[i]) for i in range(5)]

    #test vs. train seqs 1
    graph_groups([train_1, test_1], ["Train","Test"], "Train vs. Test Order 1", "", "percent accuracy", "TrainTest1.png")

    #test vs. train seqs 3
    graph_groups([train_3, test_3], ["Train","Test"], "Train vs. Test Order 3", "", "percent accuracy", "TrainTest3.png")

    #test vs. train seqs 6
    graph_groups([train_6, test_6], ["Train","Test"], "Train vs. Test Order 6", "", "percent accuracy", "TrainTest6.png")

    #order comparison
    graph_groups([test_1, test_3, test_6], ["Markov Order 1","Markov Order 3", "Markov Order 6"], "Interpolated Order Comparison", "", "percent accuracy", "OrderCompare.png")

    #transfer learning 1
    graph_groups([test_1, transfer_1], ["Amphiprion_ocellaris","Amphiprion_percula (Transfer Species)"], "Transfer Learning with Interpolated Markov Order 1", "", "percent accuracy", "Transfer1.png")

    #transfer learning 3
    graph_groups([test_3, transfer_3], ["Amphiprion_ocellaris","Amphiprion_percula (Transfer Species)"], "Transfer Learning with Interpolated Markov Order 3", "", "percent accuracy", "Transfer3.png")

    #transfer learning 6
    graph_groups([test_6, transfer_6], ["Amphiprion_ocellaris","Amphiprion_percula (Transfer Species)"], "Transfer Learning with Interpolated Markov Order 6", "", "percent accuracy", "Transfer6.png")

if __name__ == "__main__":
    main()
