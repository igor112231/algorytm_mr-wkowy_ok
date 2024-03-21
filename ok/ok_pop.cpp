#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <algorithm>
#include <cstdlib>  
#include <ctime>
#include <set>
#include <fstream>

using namespace std;

int globalBestDist=0;
string bestdna;

string generateRandomDNA(int length) {
    const string nucleotides = "ACGT";
    string dna;
    default_random_engine generator(random_device{}());
    uniform_int_distribution<int> distribution(0, 3); 

    for (int i = 0; i < length; ++i) {
        dna += nucleotides[distribution(generator)];
    }
    return dna;
}
vector<string> splitDNA(string DNA, int k)
{
    vector<string> fragments;
    for (int i = 0; i < DNA.size()-k+1; i++) {
            fragments.push_back(DNA.substr(i, k));
    }
    return fragments;
}

vector<string> makeNegativeErrors(const vector<string> &vect, float nprerror)
{
    vector<string> tempvect;
    tempvect.push_back(vect[0]);
    float size = vect.size() * nprerror;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(1, vect.size()-1);
    set<int> uniqueNumbers;

    while (uniqueNumbers.size() < size) {
        int randomNumber = distrib(gen);
        if (uniqueNumbers.find(randomNumber) == uniqueNumbers.end()) {
            uniqueNumbers.insert(randomNumber);
        }
    }

   for (auto i : uniqueNumbers)
    {
        tempvect.push_back(vect[i]);
    }
    return tempvect;
}

vector<vector<vector<int>>> creatrMatrix(const vector<vector<string>> &sortedVect)
{
    vector<vector<vector<int>>> threedvect(sortedVect.size(), vector<vector<int>>(sortedVect[0].size(), vector<int>(sortedVect[0].size())));
    for (int i = 0; i < sortedVect.size(); ++i) {
        for (int j = 0; j < sortedVect[i].size(); ++j) {
            for (int k = 0; k < sortedVect[i].size(); ++k) {
                threedvect[i][j][k] = 0;
            }
        }
    }
    for (int i = 0; i < sortedVect.size(); i++) {
        for (int j = 0; j < sortedVect[i].size(); j++) {
            for (int k = 0; k < sortedVect[i].size(); k++) {
                string str1 = sortedVect[i][j];
                string str2 = sortedVect[i][k];
                int size = sortedVect[i][k].size();
                bool stop = 0;

                for (int l = 0; l < size && stop == 0; l++)    
                {
                    string str11 = str1.substr(0, str1.length()-l);             //first letters 
                    string str22 = str2.substr(l);                              //last letters
                    if (str11 == str22)
                    {
                        threedvect[i][j][k] = l;
                        break;
                    }
                    else
                    {
                        threedvect[i][j][k] = size;
                    }
                }
                if (j == k)
                {
                    threedvect[i][j][k] = 0;
                }

            }
        }
    }
    return threedvect;
}

vector<vector<vector<float>>> creatvisMatrix(const vector<vector<vector<int>>>& distMatrix)
{
    vector<vector<vector<float>>> threedvect(distMatrix.size(), vector<vector<float>>(distMatrix[0].size(), vector<float>(distMatrix[0].size())));
    for (int i = 0; i < distMatrix.size(); ++i) {
        for (int j = 0; j < distMatrix[i].size(); ++j) {
            for (int k = 0; k < distMatrix[i].size(); ++k) {
                if (distMatrix[i][j][k] == 0)
                {
                    threedvect[i][j][k] = 0;
                }
                else
                {
                    threedvect[i][j][k] = 1 / static_cast<float>(distMatrix[i][j][k]);
                }
                if (j == k)
                {
                    threedvect[i][j][k] = 0;
                }

            }
        }
    }
    return threedvect;
}

int chooseCity(const vector<vector<float>>& pheromoneMatrix, const vector<vector<float>>& visMatrix, int counter)
{
    vector<float> probabilities;
    float sum=0;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    double randomValue = dist(gen);

    for (int i = 0; i < pheromoneMatrix[0].size(); i++) 
    {
        sum += pheromoneMatrix[i][counter] * visMatrix[i][counter] * visMatrix[i][counter];
    }

    for (int i = 0; i < pheromoneMatrix[0].size(); i++)
    {
        probabilities.push_back((pheromoneMatrix[i][counter] * visMatrix[i][counter] * visMatrix[i][counter]) / sum);

    }
    sum = 0;
    for (int i = 0; i < pheromoneMatrix.size(); i++) 
    {
        sum += probabilities[i];
    }
    for (int i = 0; i < pheromoneMatrix.size(); i++) 
    {
        probabilities[i] = probabilities[i]/sum;
    }
    for (int i = 1; i < pheromoneMatrix.size(); i++)
    {
        probabilities[i] += probabilities[i - 1];
    }

    if (probabilities[0] > randomValue)
    {
        return 0;
    }
    float sum2 = probabilities[0];
    for (int i = 1; i < probabilities.size(); i++)
    {
        sum2 += probabilities[1];
        if ((randomValue > probabilities[i-1] && randomValue <= probabilities[i]))
        {
            return i;
        }
    }
    return 1;
}

int levenshteinDist(const string& str1, const string& str2)
{
    int m = str1.length();
    int n = str2.length();
    vector<int> prevRow(n + 1, 0);
    vector<int> currRow(n + 1, 0);

    for (int j = 0; j <= n; j++) 
    {
        prevRow[j] = j;
    }

    for (int i = 1; i <= m; i++) 
    {
        currRow[0] = i;
        for (int j = 1; j <= n; j++) 
        {
            if (str1[i - 1] == str2[j - 1]) 
            {
                currRow[j] = prevRow[j - 1];
            }
            else {
                currRow[j] = 1 + min(currRow[j - 1],min(prevRow[j], prevRow[j - 1]));
            }
        }
        prevRow = currRow;
    }
    return currRow[n];
}

vector<vector<float>> runAnts(const vector<vector<float>>& visMatrix, int dnaLength, int ants, int iterations, const vector<vector<float>>& pheromoneMatrix, const vector<string> fragments, string dnaSequence)
{
    vector<vector<float>> bestpheromoneMatrix;
    int bestDistance;
    vector<int> bestpath;

    for (int iii = 0; iii < ants; iii++)
    {
        int newlength = 0;
        int counter = 0;
        int returned = 0;

        vector<vector<float>> temppheromoneMatrix;
        vector<int> path;
        string newDNA = fragments[0];
        for (const auto& row : pheromoneMatrix) {
            temppheromoneMatrix.push_back(row);
        }
        for (int i = 0; i < temppheromoneMatrix.size(); i++)
        {
            temppheromoneMatrix[0][i] = 0;
        }

        path.push_back(0);
        counter++;
        while (newlength <= dnaLength)
        {
            int itemp = 0;

            returned = chooseCity(temppheromoneMatrix, visMatrix, returned);

            path.push_back(returned);
            itemp = round(1 / visMatrix[returned][path[counter - 1]]);
            if (visMatrix[returned][path[counter - 1]] == 0)
            {
                itemp = 1;
            }
            newlength += itemp;
            newDNA += fragments[returned].substr(fragments[returned].length() - itemp);

            counter++;

            for (int i = 0; i < temppheromoneMatrix.size(); i++)
            {
                temppheromoneMatrix[returned][i] = 0;
            }
        }
        int distance = levenshteinDist(dnaSequence, newDNA);

        if ((iii == 0)||(distance < bestDistance))
        {
            for (auto& row : bestpheromoneMatrix) {
                row.clear(); 
            }
            bestpheromoneMatrix.clear();

            bestDistance = distance;
            globalBestDist = distance;
            bestdna=newDNA;
            
            for (const auto& row : temppheromoneMatrix) {
                bestpheromoneMatrix.push_back(row);
            }
            bestpath = path;
        }
    }

    for (auto& row : bestpheromoneMatrix) {
        row.clear(); 
    }
    bestpheromoneMatrix.clear();
    for (const auto& row : pheromoneMatrix) {
        bestpheromoneMatrix.push_back(row);
    }

    for (int i = 1; i < bestpath.size(); i++)
    {
        float ftemp = 1.0 / (float)bestDistance;
        bestpheromoneMatrix[bestpath[i]][bestpath[i - 1]] = bestpheromoneMatrix[bestpath[i]][bestpath[i - 1]] + ftemp + 1;
    }
    return bestpheromoneMatrix;
}

vector<vector<float>> evaporatePheromones(const vector<vector<float>>& pheromoneMatrix, float evcoef)
{
    vector<vector<float>> evpheromoneMatrix;
    for (const auto& row : pheromoneMatrix) {
        evpheromoneMatrix.push_back(row);
    }
    for (int i = 0; i < evpheromoneMatrix.size(); i++)
    {
        for (int j = 0; j < evpheromoneMatrix.size(); j++)
        {
            if(evpheromoneMatrix[i][j]<=1.0)
            { 
                evpheromoneMatrix[i][j] = 1.0;
            }
            else
            {
                evpheromoneMatrix[i][j] = evpheromoneMatrix[i][j] * evcoef;
            }
        }
    }
    return evpheromoneMatrix;
}

int main()
{
    int dnaLength = 100;
    int numDNA = 10;
    int k = 6;
    int ants = 30;
    int iterations = 30;
    float nprerror = 1;
    float evaporatecoef = 0.9;

    ofstream outputFile("data.csv");
    ofstream outputFiletxt("output.txt");
    if (outputFile.is_open())
    {
        cout << "open";
    }
    else
    {
        cout << "not";
    }
    for (int ijk = 0; ijk < 5; ijk++)
    {
        vector<string> dnaSequences;
        for (int i = 0; i < numDNA; i++) {
            dnaSequences.push_back(generateRandomDNA(dnaLength));
        }

        vector<vector<string>> dnaFragments;

        for (int i = 0; i < numDNA; i++) {
            string DNA = dnaSequences[i];
            dnaFragments.push_back(splitDNA(DNA, k));
        }

        /*for (const auto& sequenceFragments : dnaFragments) {
            for (const auto& fragment : sequenceFragments) {
                std::cout << fragment << std::endl;
            }
            std::cout << "----" << std::endl; // Separator between different sequences
        }*/

        vector<vector<string>> ndnaFragments;           //vector with errors
        if (nprerror != 1)
        {
            for (int i = 0; i < numDNA; i++) {
                vector<string> vect = dnaFragments[i];
                ndnaFragments.push_back(makeNegativeErrors(vect, nprerror));
            }
        }
        else
        {
            ndnaFragments = dnaFragments;
        }

        /*for (size_t i = 0; i < ndnaFragments.size(); ++i) {
            for (size_t j = 0; j < ndnaFragments[i].size(); ++j) {
                std::cout << ndnaFragments[i][j] << " ";
            }
            std::cout << std::endl;
        }*/
        cout << std::endl;

        vector<vector<string>> sortedVect;
        for (int i = 0; i < ndnaFragments.size(); i++) {
            sortedVect.resize(ndnaFragments.size());
            copy(ndnaFragments[i].begin(), ndnaFragments[i].end(), back_inserter(sortedVect[i]));
            sort(sortedVect[i].begin() + 1, sortedVect[i].end());
        }
        /*for (size_t i = 0; i < 1; ++i) {
            for (size_t j = 0; j < sortedVect[i].size(); ++j) {
                std::cout << sortedVect[i][j] << " ";
            }
            std::cout << std::endl;
        }*/
        vector<vector<vector<int>>> distMatrix = creatrMatrix(sortedVect);

        /*for (size_t i = 0; i < distMatrix.size(); ++i) {
            std::cout << "Matrix " << i + 1 << ":\n";
            for (size_t j = 0; j < distMatrix[i].size(); ++j) {
                for (size_t l = 0; l < distMatrix[i][j].size(); ++l) {
                    std::cout << distMatrix[i][j][l] << " ";
                }
                std::cout << std::endl; // New line at the end of each row
            }
            std::cout << std::endl; // Separate each matrix with a new line
        }*/

        vector<vector<vector<float>>> visMatrix = creatvisMatrix(distMatrix);
        /*for (size_t i = 0; i < visMatrix.size(); ++i) {
            std::cout << "Matrix " << i + 1 << ":\n";
            for (size_t j = 0; j < visMatrix[i].size(); ++j) {
                for (size_t l = 0; l < visMatrix[i][j].size(); ++l) {
                    std::cout << visMatrix[i][j][l] << " ";
                }
                std::cout << std::endl; // New line at the end of each row
            }
            std::cout << std::endl; // Separate each matrix with a new line
        }*/

        vector<int> allLewDist;
        vector<string>allDNA;

        for (int iii = 0; iii < numDNA; iii++)
        {
            vector<vector<float>> pheromoneMatrix(visMatrix[0].size(), vector <float>(visMatrix[0].size(), 1));

            for (int i = 0; i < iterations; i++)
            {
                pheromoneMatrix = runAnts(visMatrix[iii], dnaLength, ants, iterations, pheromoneMatrix, sortedVect[iii], dnaSequences[iii]);
                pheromoneMatrix = evaporatePheromones(pheromoneMatrix, evaporatecoef);
            }
            allLewDist.push_back(globalBestDist);
            allDNA.push_back(bestdna);
        }
        float mean = 0.0;
        for (int i = 0; i < allLewDist.size(); i++)
        {
            mean += allLewDist[i];
            cout << allLewDist[i] << " ";
        }
        mean = mean / allLewDist.size();
        cout << endl << mean;
        outputFile << mean << endl;
        

        outputFiletxt << "evaporation coeffincient = " << evaporatecoef << endl;
        for (int i = 0; i < numDNA; i++) {
            outputFiletxt << "original DNA: " << dnaSequences[i] << endl << "new DNA: " << allDNA[i] << endl << "Levenstein distance: " << allLewDist[i] << endl;
        }
        evaporatecoef -= 0.02;

    }
}
