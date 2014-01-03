EDF-Reader
==========

Allows users to read .edf EEG files in Matlab.

Reading small EDF files in Matlab is relatively easy and there are multiple methods available online to do so. However, when the EDF files are very large, these tools do not work nicely as they usually try to read the entire file into memory.

This EDF-Reader creates an object representing the EDF file, and provides random access to the data in the EDF-file. When you create the object, the EDF header is read into the object and the GETDATA method can be used to read (part) of the available data.

    Example: Reading the first 10000 values of the second channel in 'newFile.edf'
      >> out = EDFData('newFile.edf')
        out = 
          EDFData with properties:
                  ver: 0
            patientID: [1x80 char]
             recordID: [1x80 char]
            startDate: [1 1 70]
            startTime: [0 0 0]
                bytes: 4352
              records: 60
             duration: 0.9982
                   ns: 16
                label: {1x16 cell}
           transducer: {1x16 cell}
                units: {1x16 cell}
            prefilter: {1x16 cell}
              samples: [1x16 double]
                   sf: [1x16 double]
           conversion: [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
             dcOffset: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
             
      >> data = out.getData([1 10000],2);
          data =
               -23
                59
               111
               ...
               
**Install:** 
Add the @EDFData folder to the Matlab path.
