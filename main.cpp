#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

//namespace variables
using std::ifstream;
using std::cout;
using std::cin;
using std::string;

int main(int argc, char *argv[])
{
	if (argc < 5)
	  return 1;
		
	//Check arguments that are passed in
	cout << argc << "\n";
	cout << (argc - 1) / 2 << "\n";
	
	//Amount of runs to check; each run has two files. First argument is program name
	int fileCount = (argc - 1) / 2;
	
	//Define array variables
	int pointCount, fpointCount;
	double ***positions, ***positionsNext;
	int **ids, **idsNext;
	int **remapPoints;
	int **remapPointsNext, **mapPointsNext;
	
	//Initialize array variables
	positions = new double**[fileCount];
	positionsNext = new double**[fileCount];
	remapPoints = new int*[fileCount];
	mapPointsNext = remapPointsNext = new int*[fileCount];
	ids = new int*[fileCount];
	idsNext = new int*[fileCount];
	
	//Initialize array of filestreams
	ifstream *myFiles = new ifstream[fileCount];
	ifstream *myFilesNext = new ifstream[fileCount];
	
	{ //First pair of files
	  int i = 0;
		cout << "---First File!---\n"; //Output notification of where we are
    string testStr, testStr2; //Define two strings
	  myFiles[i].open(argv[(2*i)+1]);     //Initialize filestream to first file in a pair
	  myFilesNext[i].open(argv[(2*i)+2]); //Initialize filestream to second file in a pair
		
		
		for (testStr = ""; testStr != "Piece";  )
		{ //Loop over file, until next word is "Piece"
			myFiles[i].ignore(1024, '<'); //Skip <
			myFiles[i] >> testStr; //Read next word
		}
		for (testStr = ""; testStr != "Piece";  )
		{ //Same as above, for second file
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		
		
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{ //Move to end of 'NumberOfPoints="'
			myFiles[i].ignore(1024, ' ');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFilesNext[i].ignore(1024, ' ');
			myFilesNext[i] >> testStr;
		}
		
		
		fpointCount = stoi(testStr.substr(16,16),nullptr);
		//cout << stoi(testStr.substr(16,16),nullptr) << "\n";
		
		//Initialize arrays for positions and IDs, since we now know the locations
		positions[i] = new double*[fpointCount];     //Initial positions
		positionsNext[i] = new double*[fpointCount]; //Ending positions
		remapPoints[i] = new int[fpointCount];
		mapPointsNext[i] = new int[fpointCount];
		ids[i] = new int[fpointCount];
		idsNext[i] = new int[fpointCount];
		
		
		for (testStr = ""; testStr != "DataArray";  )
		{ //Skip ahead until we get to "DataArray"
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr != "DataArray";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		
		for (int j = 0; j < fpointCount; j++)
		{ 
		  //Remappoints[i][j] = k:
			//means jth particle in the ith file, corresponds to the kth particles in the 1st file
			//since this is the first file, j = k
		  remapPoints[i][j] = j;
			
			myFiles[i].ignore(1024, '\n'); //Skips a line to data, or go to next line
			myFilesNext[i].ignore(1024, '\n');
		  positions[i][j] = new double[3];
		  positionsNext[i][j] = new double[3];
			myFiles[i] >> positions[i][j][0];  //Read in position along each axis
			myFiles[i] >> positions[i][j][1];
			myFiles[i] >> positions[i][j][2];
			myFilesNext[i] >> positionsNext[i][j][0]; //Same for end file
			myFilesNext[i] >> positionsNext[i][j][1];
			myFilesNext[i] >> positionsNext[i][j][2];
			
			//Output positions to user
			cout << " Position [" << i << "][" << j << 
			  "]: " << positions[i][j][0] <<
			  ", "  << positions[i][j][1] << 
			  ", "  << positions[i][j][2] << "\n";
			cout << " Next Pos [" << i << "][" << j << 
			  "]: " << positionsNext[i][j][0] <<
			  ", "  << positionsNext[i][j][1] << 
			  ", "  << positionsNext[i][j][2] << "\n";
		}
		
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{ //Skip to ID section
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFiles[i].ignore(1024, '<');
				myFiles[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFiles[i].ignore(1024, ' ');
				myFiles[i] >> testStr;
			}
			
		}
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{ //outer loop checks that name is indeed "ID"
			for (testStr = ""; testStr != "DataArray";  )
			{ //Inner loop skips to dataArray section
				myFilesNext[i].ignore(1024, '<');
				myFilesNext[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{ //Inner loop skips to Name section
				myFilesNext[i].ignore(1024, ' ');
				myFilesNext[i] >> testStr;
			}
		}
		
		
		for (int j = 0; j < fpointCount; j++)
		{ //Input IDs 
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');
			myFiles[i] >> ids[i][j];
			myFilesNext[i] >> idsNext[i][j];
			cout << " id [" << i << "][" << j << 
			  "]: " << ids[i][j] <<
			  ", "  << idsNext[i][j] << "\n";
		}
		
		
		for (int j = 0; j < fpointCount; j++)
		{ //If the order of the "next" file is different from the first file, re-order according to IDs
		
			//mapPointsNext[i][j] = k
			//In ith file pair, jth particle in the first file corresponds to the kth particle in the second file
		
		  if (ids[i][j] == idsNext[i][j])
			  mapPointsNext[i][j] = j;
			else
			{
			  for (int k = 0; k < fpointCount; k++)
				{
				  if (ids[i][j] == idsNext[i][k])
					  mapPointsNext[i][j] = k;
				}
			}
		}
		
		cout << "\n";
	}
	
	
	//Repeat same as above for the other files. fpointCount is unique to first file,
	//but the pointCount of subsequent files doesn't matter
	for (int i = 1; i < fileCount; i++)
	{
		cout << "----File [" << i << "]!----\n";
	  string testStr;
	  myFiles[i].open(argv[(2*i)+1]);
	  myFilesNext[i].open(argv[(2*i)+2]);
		for (testStr = ""; testStr != "Piece";  )
		{
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr != "Piece";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFiles[i].ignore(1024, ' ');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFilesNext[i].ignore(1024, ' ');
			myFilesNext[i] >> testStr;
		}
		
		// cout << "Test1\n";
		
		pointCount = stoi(testStr.substr(16,16),nullptr);
		//cout << stoi(testStr.substr(16,16),nullptr) << "\n";
		positions[i] = new double*[pointCount];
		positionsNext[i] = new double*[pointCount];
		remapPoints[i] = new int[fpointCount];
		mapPointsNext[i] = new int[pointCount];
		ids[i] = new int[pointCount];
		idsNext[i] = new int[pointCount];
		
		for (int j = pointCount; j < fpointCount; j++)
		{
		  remapPoints[i][j] = -1;
		}
		
		for (testStr = ""; testStr != "DataArray";  )
		{
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr != "DataArray";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		// cout << "test2\n";
		for (int j = 0; j < pointCount; j++)
		{
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');
		  positions[i][j] = new double[3];
		  positionsNext[i][j] = new double[3];
			myFiles[i] >> positions[i][j][0];
			myFiles[i] >> positions[i][j][1];
			myFiles[i] >> positions[i][j][2];
			myFilesNext[i] >> positionsNext[i][j][0];
			myFilesNext[i] >> positionsNext[i][j][1];
			myFilesNext[i] >> positionsNext[i][j][2];
			
			
			///Difference here///
			///               ///
			///_______________///
			
			if (j < fpointCount)
			{	//For every jth particle in the ith file
				if (positions[i][j][0] == positions[0][j][0] &&
						positions[i][j][1] == positions[0][j][1] &&
						positions[i][j][2] == positions[0][j][2])
				{ //Check if position is the same as position of jth particle in first file
					remapPoints[i][j] = j;
				}
				else
				{
					for (int k = 0; k < fpointCount; k++)
					{ //For every kth particle in the first
						if (positions[i][j][0] == positions[0][k][0] &&
								positions[i][j][1] == positions[0][k][1] &&
								positions[i][j][2] == positions[0][k][2])
						{ //When we find the kth particle that is in the same position as the jth particle in the first file, save that info
							if (k < fpointCount)
								remapPoints[i][k] = j;
						}
					}
				}
			}
			//cout << "test2.1\n";
		}
		
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFiles[i].ignore(1024, '<');
				myFiles[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFiles[i].ignore(1024, ' ');
				myFiles[i] >> testStr;
			}
			
		}
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFilesNext[i].ignore(1024, '<');
				myFilesNext[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFilesNext[i].ignore(1024, ' ');
				myFilesNext[i] >> testStr;
			}
		}
		
		for (int j = 0; j < pointCount; j++)
		{
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');
			myFiles[i] >> ids[i][j];
			myFilesNext[i] >> idsNext[i][j];
			cout << " id [" << i << "][" << j << 
			  "]: " << ids[i][j] <<
			  ", "  << idsNext[i][j] << "\n";
		}
		
		
		for (int j = 0; j < pointCount; j++)
		{
		  if (ids[i][j] == idsNext[i][j])
			  mapPointsNext[i][j] = j;
			else
			{
			  for (int k = 0; k < fpointCount; k++)
				{
				  if (ids[i][j] == idsNext[i][k])
					  mapPointsNext[i][j] = k;
				}
			}
		}
		
		
		for (int j = 0; (j < pointCount) && (j < fpointCount); j++)
		{
			cout << " Position [" << i << "][" << (remapPoints[i][j]) << 
			  "]: " << positions[i][(remapPoints[i][j])][0] <<
			  ", "  << positions[i][(remapPoints[i][j])][1] << 
			  ", "  << positions[i][(remapPoints[i][j])][2] <<
        ". Remap[" << i << "][" << j << "]: " << remapPoints[i][j] << "\n";
			cout << " Next Pos [" << i << "][" << (remapPoints[i][j]) << 
			  "]: " << positionsNext[i][(remapPoints[i][j])][0] <<
			  ", "  << positionsNext[i][(remapPoints[i][j])][1] << 
			  ", "  << positionsNext[i][(remapPoints[i][j])][2] <<
        ". Remap[" << i << "][" << j << "]: " << remapPoints[i][j] << "\n";
		}
		cout << "\n";
	}
	
	
	//Comparisons start here
	cout << "\n---Compare Begin---\n\n";
	
	for (int i = 0; i < fpointCount; i++)
	{
	  /*cout << "file[0]:Point[" << i << "]: <" << 
      positionsNext[0][mapPointsNext[0][i]][0] << "," <<
      positionsNext[0][mapPointsNext[0][i]][1] << "," <<
      positionsNext[0][mapPointsNext[0][i]][2] << "," <<
		  ">; ID: " << idsNext[0][mapPointsNext[0][i]] <<"\n";*/
			cout << "Next Point\n";
			int k2 = mapPointsNext[0][i]; //index of End position in file J-1 of ith particle
			double errorOld = 0;
	  for (int j = 1; j < fileCount; j++)
		{
		  int k1 = remapPoints[j][i]; //index of Start position in file J for ith particle
			k1 = mapPointsNext[j][k1]; //Now index of end position in file J for ith particle
		  //int k2 = remapPoints[j-1][i];
			//k2 = mapPointsNext[j-1][k1];
			double errX = positionsNext[j][k1][0] - positionsNext[j-1][k2][0],
			       errY = positionsNext[j][k1][1] - positionsNext[j-1][k2][1],
			       errZ = positionsNext[j][k1][2] - positionsNext[j-1][k2][2];
	    double error = sqrt((errX * errX) + (errY * errY) + (errZ * errZ));
			
		  /*cout << " file[" << j-1 << "].point[" << k1 << "]: <" << 
				positionsNext[j-1][k2][0] << "," <<
				positionsNext[j-1][k2][1] << "," <<
				positionsNext[j-1][k2][2] << "," <<
				">; ID: " << idsNext[j-1][k2] << "\n";*/
			
			double errorLg = log2(error/errorOld);
			cout << "                  Order: {" << errorLg << "}\n";
			
			cout << " [" << j-1 << "],[" << j << "]\n";
			 // ID: (" << idsNext[j-1][k2]
			     // << "," << idsNext[j][k1] << ")\n";
			cout << "  error: " << error << "\n";
			
			errorOld = error;
			
		  /*cout << " file[" << j << "].point[" << k1 << "]: <" << 
				positionsNext[j][k1][0] << "," <<
				positionsNext[j][k1][1] << "," <<
				positionsNext[j][k1][2] << "," <<
				">; ID: " << idsNext[j][k1] << "\n";*/
			k2 = k1;
			//cout << "\n";
		}
		cout << "\n";
	}
	
	return 0;
}
