classdef xnodes < factor_node
% EVIDENT_NODE derived from factor_node class provides a message to other
% nodes corresponding an observed event.
%
% Description of methods:
%   EVIDENT_NODE    instantiate an evident node
%   SETUP_INIT_MSG  send out initial message
%   RESET           reset evident node
%   SETUP_LINK      setup linklist
%
% Usage:
% Initialize an evident nodes m and a factor node n, with id = 1,2 respectively
%     m = evident_node(1);
%     n = factor_node(2);
% Set up the linklist such that m,n,o are linked
%     m.setup_link({n});
%     n.setup_link({m});
% Set up the initial message in row vector.  Each element represents a
% probability.
% For binary node, the 1st/2nd element are event=false and event=true
%     m.setup_init_msg([0.9 0.1]);
% Reset nodes
%     m.reset();
%     n.reset();

properties
       init_msg
       known
   end

   % The evident node sends out a message upon
   % receiving the initial message.
   % It takes no action on messages from other nodes.
   methods
       function s = xnodes(id)
           s = s@factor_node(id);
           s.init_msg = {};
           s.known = 0;
       end

       function setup_init_msg(s, msg)
           s.known = 1;
           if size(msg,1) ~= 1
               sprintf('Node %0.0f ERROR from setup_init_msg: message must be row vector !!!\n',s.id)
           else
               s.init_msg{1} = msg;
               % construct message
               snd_msg = s.factor_fun(s.init_msg,{},{});
               % send message
               for i=1:length(s.linklist)
                   s.linklist{i}.rx_msg(s, snd_msg);
                   % update outbound record
                   s.to_node{i} = s.linklist{i};
                   s.outbound_msg{i} = snd_msg;
                   s.to_id = [s.to_id s.linklist{i}.id];
               end
           end
       end
       function reset(s)
           reset@factor_node(s);
           s.init_msg = {};
       end
       function setup_link(s,linklist)
           setup_link@factor_node(s,linklist);
       end
%    end
% 
%    methods
       function msg = factor_fun(s,in_msg,from_id,to_id)
           in_msg{1}
           if length(in_msg) == 1
               msg = in_msg{1};
           else
               msg = [1 1 1];
               for i=1:length(in_msg)
                   msg = msg .* in_msg{i};
               end
               % re-scale msg
               sum_msg = sum(msg);
               if sum_msg == 0
                   msg = ones(size(msg));
               end
               msg = msg/sum(msg);
           end
       end
   end
       
   methods (Access = protected)
       function rx_msg(s, from_node, msg)
           s.from_node{1} = from_node;
           s.inbound_msg{1} = msg;
           s.from_id = [from_node.id];
       end
   end
end 
