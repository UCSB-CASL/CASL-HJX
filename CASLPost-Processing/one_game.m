function winnings = one_game(num_other_players)
%ONE_GAME  Simulate a single Massachusetts WinFall drawing for Jerry
%
%   winnings = ONE_GAME(K) returns Jerry’s net winnings (in dollars) when
%   he buys 200 000 tickets and there are K *other* players, each of whom
%   also buys 200 000 tickets.
%
%   Payout rules follow the project statement exactly:
%     • Jackpot = 2 000 000 + 240 000 × K
%     • If ≥1 ticket matches 6 numbers  →  jackpot split evenly.
%     • Otherwise the jackpot rolls down:
%           – add 26 % to every 5-match ticket
%           – add 47 % to every 4-match ticket
%           – add 27 % to every 3-match ticket
%     • Base payouts: 5-match   $4 000
%                     4-match   $150
%                     3-match   $5
%                     2-match   $2   (value of free ticket)
%     • Jerry’s net = total payout − $400 000 (cost of his tickets).
%
%   Requires NUM_MATCH.M and ONE_PLAYER.M.

% -------------------------------------------------------------------------
% 1.  Set up the drawing
% -------------------------------------------------------------------------
jackpot = 2000000 + 240000 * num_other_players;   % total jackpot (USD)
winner  = randperm(46,6);                            % the six winning balls

% Jerry’s 200 000 tickets
[m2J,m3J,m4J,m5J,m6J] = one_player(winner);

% Other players (aggregate counts)
m2O=0; m3O=0; m4O=0; m5O=0; m6O=0;
for p = 1:num_other_players
    [a2,a3,a4,a5,a6] = one_player(winner);
    m2O = m2O + a2;  m3O = m3O + a3;  m4O = m4O + a4;
    m5O = m5O + a5;  m6O = m6O + a6;
end

% Totals across *all* players (needed for split/roll-down logic)
tot2 = m2J + m2O;
tot3 = m3J + m3O;
tot4 = m4J + m4O;
tot5 = m5J + m5O;
tot6 = m6J + m6O;

% -------------------------------------------------------------------------
% 2.  Compute Jerry’s gross payout
% -------------------------------------------------------------------------
payout = 0;

if tot6 > 0
    % ------------------------------------------------- Jackpot is hit
    jackpot_each = jackpot / tot6;
    payout = payout + m6J * jackpot_each;          % Jerry’s share
    % standard payouts for lower tiers
    payout = payout + 4000 * m5J;
    payout = payout +   150 * m4J;
    payout = payout +     5 * m3J;
else
    % ------------------------------------------------- Roll-down
    % 5-match tickets
    if tot5 > 0
        payout_each5 = 4000 + 0.26 * jackpot / tot5;
        payout = payout + m5J * payout_each5;
    end
    % 4-match tickets
    if tot4 > 0
        payout_each4 =   150 + 0.47 * jackpot / tot4;
        payout = payout + m4J * payout_each4;
    end
    % 3-match tickets
    if tot3 > 0
        payout_each3 = 5 + 0.27 * jackpot / tot3;
        payout = payout + m3J * payout_each3;
    end
end

% 2-match tickets (no roll-down component, always $2)
payout = payout + 2 * m2J;

% -------------------------------------------------------------------------
% 3.  Subtract Jerry’s ticket outlay
% -------------------------------------------------------------------------
cost     = 200000 * 2;     % 200 k tickets at $2 each
winnings = payout - cost;
end
