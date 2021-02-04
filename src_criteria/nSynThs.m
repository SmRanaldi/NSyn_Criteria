function [nSyn]=nSynThs(VAFIn,ths)

nSyn=min(find(VAFIn>ths));