%% Implicit expansion
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
% 
%All parameters and results are in SI base units.
%
%
%
%This class contains some auxiliary functions to facilitate the use of
%implicit expansion.
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14


classdef implExp
    methods(Static)
        %Predict the size of an array when using implicit expansion
        function sz=size(varargin)
            dims=max(cellfun(@(x) ndims(x),varargin));
            szCell=cellfun(@(x) implExp.padsz(dims,x),varargin,'UniformOutput',false);
            szs=cell2mat(szCell');

            %all dimensions not equal to 1 must be the same
            not1=szs~=1;
            check=arrayfun(@(x) all(szs(not1(:,x),x)==szs(find(not1(:,x),1),x)),1:size(szs,2));

            if all(check)
                sz=max(szs,[],1);
                sz(any(szs==0,1))=0;
            else
                throw(MException('implExp:incompatibleArrays','Arrays have incompatible sizes'));
            end
        end


        %Same as "size", but with one dimension discretized
        function sz=sizeDisc(reg,disc,discDim)
            %reg        regular variables, not discretized
            %disc       discretized variables
            %discDim    dimension in which the disc-variables are discretized

            discSz=implExp.szArr(disc{:});
            discSz(:,discDim)=1;

            param=arrayfun(@(x) NaN([discSz(x,:)]),1:length(disc),'UniformOutput',false);

            sz=implExp.size(reg{:},param{:});
        end


        %Transform all variable into row vectors of equal length
        function varargout=normalize(sz,varargin)
            n=prod(sz);
            if n~=0
                pad=@(x) implExp.padsz(numel(sz),x);
                varargout=cellfun(@(x) reshape(repmat(x,sz./pad(x)),1,n),varargin,'UniformOutput',false);
            else
                varargout=repmat({NaN(sz)},1,length(varargin));
            end
        end


        %Same as "normalize", but with one dimension discretized
        function [reg,disc]=normDisc(sz,reg,disc,discDim)
            %reg        regular variables, not discretized
            %disc       discretized variables
            %discDim    dimension in which the disc-variables are discretized

            disc=cellfun(@(x) shiftdim(x,discDim-1),disc,'UniformOutput',false);

            c=[disc,reg];
            expand=implExp.szArr(c{:});
            expand(1:length(disc),1)=1;

            disc=arrayfun(@(x) ...
                    reshape(repmat(disc{x},sz./expand(x,:)),size(disc{x},1),prod(sz)),...
                    1:length(disc),'UniformOutput',false);

            [reg{:}]=implExp.normalize(sz,reg{:});
        end
    end
    
    
    methods(Static, Access=private)
        function sz=szArr(varargin)
            dims=max(cellfun(@(x) ndims(x),varargin));
            sz=cell2mat(cellfun(@(x) implExp.padsz(dims,x),varargin,'UniformOutput',false)');
        end


        function padsz=padsz(dims,x)
            padsz=[size(x),ones(1,dims-ndims(x))];
        end
    end
end




