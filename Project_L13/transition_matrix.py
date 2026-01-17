import json
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path

class TransitionMatrixGenerator:
    """Generate transition matrices from sequences and save as JSON"""
    
    @staticmethod
    def dna_transition_matrix(sequence: str) -> dict:
        """
        Compute transition probabilities between DNA letters (A, T, G, C).
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            Dictionary with transition matrix and metadata
        """
        # Clean sequence - keep only ATGC
        sequence = sequence.upper()
        sequence = ''.join(c for c in sequence if c in 'ATGC')
        
        if len(sequence) < 2:
            raise ValueError("Sequence too short for transition analysis")
        
        letters = ['A', 'T', 'G', 'C']
        transitions = defaultdict(lambda: defaultdict(int))
        
        # Count transitions
        for i in range(len(sequence) - 1):
            current = sequence[i]
            next_letter = sequence[i + 1]
            transitions[current][next_letter] += 1
        
        # Convert to probabilities
        matrix = {}
        for letter in letters:
            matrix[letter] = {}
            total = sum(transitions[letter].values())
            
            for next_letter in letters:
                if total > 0:
                    matrix[letter][next_letter] = transitions[letter][next_letter] / total
                else:
                    matrix[letter][next_letter] = 0.0
        
        return {
            "type": "DNA Sequence Transitions",
            "sequence_length": len(sequence),
            "sequence_sample": sequence[:50] + "..." if len(sequence) > 50 else sequence,
            "letters": letters,
            "transition_matrix": matrix,
            "description": "Rows: current letter, Columns: next letter probability"
        }
    
    @staticmethod
    def text_transition_matrix(text: str) -> dict:
        """
        Compute transition probabilities between words in English text.
        
        Args:
            text: English text string
        
        Returns:
            Dictionary with transition matrix and metadata
        """
        # Split into words and clean
        words = text.split()
        words = [w.lower() for w in words if w]
        
        if len(words) < 2:
            raise ValueError("Text too short for transition analysis")
        
        transitions = defaultdict(lambda: defaultdict(int))
        word_list = list(set(words))  # Unique words
        
        # Count transitions between words
        for i in range(len(words) - 1):
            current_word = words[i]
            next_word = words[i + 1]
            transitions[current_word][next_word] += 1
        
        # Convert to probabilities
        matrix = {}
        for word in word_list:
            matrix[word] = {}
            total = sum(transitions[word].values())
            
            for next_word in word_list:
                if total > 0:
                    matrix[word][next_word] = transitions[word][next_word] / total
                else:
                    matrix[word][next_word] = 0.0
        
        return {
            "type": "Word Sequence Transitions",
            "text_length": len(text),
            "total_words": len(words),
            "unique_words": len(word_list),
            "word_list": word_list,
            "text_sample": text[:100] + "..." if len(text) > 100 else text,
            "transition_matrix": matrix,
            "description": "Rows: current word, Columns: next word probability"
        }
    
    @staticmethod
    def save_to_json(data: dict, filename: str) -> None:
        """Save transition matrix data to JSON file"""
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"✓ Saved to {filename}")


def main():
    # ==================== TASK 1: DNA SEQUENCE ====================
    print("=" * 70)
    print("TASK 1: DNA SEQUENCE TRANSITION MATRIX")
    print("=" * 70)
    
    # Random DNA sequence (about 50 letters)
    dna_sequence = "ATGCATGCTAGCTAGCTAGCTAATGCTAGCTAGCTAGCTAGCTAATGCTAGCTA"
    print(f"\nDNA Sequence ({len(dna_sequence)} letters):")
    print(dna_sequence)
    
    # Generate transition matrix
    dna_result = TransitionMatrixGenerator.dna_transition_matrix(dna_sequence)
    
    print("\nTransition Matrix (DNA):")
    for letter, transitions in dna_result["transition_matrix"].items():
        print(f"\n{letter} ->", end=" ")
        for next_letter, prob in transitions.items():
            print(f"{next_letter}: {prob:.3f}", end="  ")
    
    # Save to JSON
    TransitionMatrixGenerator.save_to_json(dna_result, "dna_transition_matrix.json")
    
    
    # ==================== TASK 2: ENGLISH TEXT ====================
    print("\n\n" + "=" * 70)
    print("TASK 2: WORD SEQUENCE TRANSITION MATRIX")
    print("=" * 70)
    
    # Random English text (about 300 characters)
    english_text = """The quick brown fox jumps over the lazy dog. The dog was sleeping
    under the big tree. The tree was very old and tall. The fox and the dog became
    friends after that day. They would play together every morning. The morning was
    always beautiful and peaceful. The peaceful environment made them happy. They
    were happy and content with their new friendship."""
    
    print(f"\nEnglish Text ({len(english_text)} characters):")
    print(english_text[:150] + "...")
    
    # Generate transition matrix
    text_result = TransitionMatrixGenerator.text_transition_matrix(english_text)
    
    print(f"\nTotal words: {text_result['total_words']}")
    print(f"Unique words: {text_result['unique_words']}")
    
    print("\nWord List:")
    print(", ".join(text_result["word_list"][:10]), "...")
    
    print("\nTransition Matrix (Sample - first 3 words):")
    for word in list(text_result["transition_matrix"].keys())[:3]:
        print(f"\n{word} ->", end=" ")
        for next_word, prob in list(text_result["transition_matrix"][word].items())[:3]:
            if prob > 0:
                print(f"{next_word}: {prob:.3f}", end="  ")
    
    # Save to JSON
    TransitionMatrixGenerator.save_to_json(text_result, "text_transition_matrix.json")
    
    
    # ==================== CUSTOM INPUT OPTION ====================
    print("\n\n" + "=" * 70)
    print("CUSTOM INPUT OPTION")
    print("=" * 70)
    
    choice = input("\nWould you like to analyze your own data? (y/n): ").lower()
    
    if choice == 'y':
        data_type = input("DNA sequence or Text? (dna/text): ").lower()
        
        if data_type == 'dna':
            seq = input("Enter your DNA sequence (ATGC only): ")
            try:
                result = TransitionMatrixGenerator.dna_transition_matrix(seq)
                TransitionMatrixGenerator.save_to_json(result, "custom_dna_matrix.json")
                print("\nTransition Matrix:")
                for letter, trans in result["transition_matrix"].items():
                    print(f"{letter}: {trans}")
            except ValueError as e:
                print(f"Error: {e}")
        
        elif data_type == 'text':
            text = input("Enter your English text: ")
            try:
                result = TransitionMatrixGenerator.text_transition_matrix(text)
                TransitionMatrixGenerator.save_to_json(result, "custom_text_matrix.json")
                print(f"\nWords found: {result['total_words']}")
                print("Transition Matrix saved!")
            except ValueError as e:
                print(f"Error: {e}")


if __name__ == "__main__":
    main()