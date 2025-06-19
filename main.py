import numpy as np

def create_score_matrix(seq1, seq2, match=3, mismatch=-3, gap_penalty=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)
    
    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                diagonal = score_matrix[i-1][j-1] + match
            else:
                diagonal = score_matrix[i-1][j-1] + mismatch
            
            up = score_matrix[i-1][j] + gap_penalty
            left = score_matrix[i][j-1] + gap_penalty
            
            score = max(0, diagonal, up, left)
            score_matrix[i][j] = score
            
            if score > max_score:
                max_score = score
                max_pos = (i, j)
    
    return score_matrix, max_pos

def traceback(score_matrix, seq1, seq2, max_pos, match=3, mismatch=-3, gap_penalty=-2):
    i, j = max_pos
    align1 = []
    align2 = []
    
    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        diagonal = score_matrix[i-1][j-1]
        up = score_matrix[i-1][j]
        left = score_matrix[i][j-1]
        
        if current_score == diagonal + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif current_score == up + gap_penalty:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        elif current_score == left + gap_penalty:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    return ''.join(reversed(align1)), ''.join(reversed(align2))

def print_alignment(align1, align2):
    print("\nBest Local Alignment:")
    print("Sequence 1:", align1)
    print("Sequence 2:", align2)
    
    match_line = []
    for a, b in zip(align1, align2):
        if a == b:
            match_line.append('|')
        else:
            match_line.append(' ')
    print("          ", ''.join(match_line))

def load_database(filename="database.txt"):
    try:
        with open(filename, 'r') as file:
            return [line.strip().upper() for line in file if line.strip()]
    except FileNotFoundError:
        print(f"{filename} file not found. Using sample database.")
        return [
            "ATCGATCGATCG",
            "GGGATCCCCGAT",
            "TTTAGCTAGCTA",
            "CATCGAATCGGG",
            "AAATTTCCGGGA"
        ]

def save_database(database, filename="database.txt"):
    with open(filename, 'w') as file:
        for seq in database:
            file.write(seq + '\n')

def main():
    print("Mini BLAST Tool - Local Alignment (Smith-Waterman)")
    
    database = load_database()
    
    while True:
        print("\nMenu:")
        print("1. Align sequence")
        print("2. Show database")
        print("3. Add sequence to database")
        print("4. Exit")
        
        choice = input("Your choice (1-4): ")
        
        if choice == "1":
            print("\nSequences in Database:")
            for i, seq in enumerate(database, 1):
                print(f"{i}. {seq}")
            
            user_seq = input("\nEnter your sequence (A,T,C,G only): ").upper()
            
            valid_chars = {'A', 'T', 'C', 'G'}
            if not all(c in valid_chars for c in user_seq):
                print("Error: Sequence must contain only A, T, C, G characters.")
                continue
            
            best_score = -1
            best_alignment = ("", "")
            best_db_seq = ""
            
            for db_seq in database:
                score_matrix, max_pos = create_score_matrix(user_seq, db_seq)
                align1, align2 = traceback(score_matrix, user_seq, db_seq, max_pos)
                
                current_score = score_matrix[max_pos[0]][max_pos[1]]
                if current_score > best_score:
                    best_score = current_score
                    best_alignment = (align1, align2)
                    best_db_seq = db_seq
                    best_matrix = score_matrix
                    best_max_pos = max_pos
            
            print("\nHighest Score:", best_score)
            print("Database Sequence:", best_db_seq)
            print("Your Sequence:", user_seq)
            
            print("\nScore Matrix:")
            print(best_matrix)
            
            print_alignment(best_alignment[0], best_alignment[1])
        
        elif choice == "2":
            print("\nDatabase Contents:")
            for i, seq in enumerate(database, 1):
                print(f"{i}. {seq}")
        
        elif choice == "3":
            new_seq = input("\nEnter new sequence to add (A,T,C,G only): ").upper()
            
            valid_chars = {'A', 'T', 'C', 'G'}
            if not all(c in valid_chars for c in new_seq):
                print("Error: Sequence must contain only A, T, C, G characters.")
                continue
            
            database.append(new_seq)
            save_database(database)
            print(f"'{new_seq}' added to database.")
        
        elif choice == "4":
            # Exit program
            print("Exiting program...")
            break
        
        else:
            print("Invalid choice. Please enter a number between 1-4.")

if __name__ == "__main__":
    main()
