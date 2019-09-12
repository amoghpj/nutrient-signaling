import yaml

class TimeCourse:
    self.data = []
    
    def read_data(self, path_to_yaml):
        with open(path_to_yaml,'r') as infile:
            self.data = yaml.safe_load(infile)
