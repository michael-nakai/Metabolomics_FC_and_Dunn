import pandas as pd
from helper import ImportedData


def main():

    # Define filepaths to caecum_all and liver_all, along with the outputs
    caecum_path = r'C:\Users\Michael\Desktop\Code\Python\Rikeish_Metabolomics\raw-data\Caecum_ALL.csv'
    liver_path = r'C:\Users\Michael\Desktop\Code\Python\Rikeish_Metabolomics\raw-data\Liver_ALL.csv'
    caecum_output = r'C:\Users\Michael\Desktop\Code\Python\Rikeish_Metabolomics\outputs\Caecum_FC.csv'
    liver_output = r'C:\Users\Michael\Desktop\Code\Python\Rikeish_Metabolomics\outputs\Liver_FC.csv'

    # Read the main csv files
    print("Reading files")
    caecum_all_pd = pd.read_csv(caecum_path)
    liver_all_pd = pd.read_csv(liver_path)

    # Construct two ImportedData objects
    print("Constructing objects")
    caecum = ImportedData(caecum_all_pd)
    liver = ImportedData(liver_all_pd)

    # Calculate values
    print("Calculating caecum values")
    caecum.calculate_values()
    print("Calculating liver values")
    liver.calculate_values()

    # Save results
    print("Saving results")
    caecum.save_results(caecum_output)
    liver.save_results(liver_output)


if __name__ == '__main__':
    main()
