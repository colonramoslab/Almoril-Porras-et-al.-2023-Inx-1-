function sc = fromString(str, varargin)
%converts from a string representation to the array of segmentation
%clusters
%function sc = fromString(str, eset)
import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
obj = yaml.load(str);
sc = SegmentationCluster.fromHashMap(obj, varargin{:});
